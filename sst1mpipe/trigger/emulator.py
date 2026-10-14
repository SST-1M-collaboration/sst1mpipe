"""Software emulation of the SST-1M DigiCam trigger, from R0 data.

Two emulators share the base class ``TriggerEmulator``:

* ``TriggerPatch7Emulator``, the trigger running on the telescopes: sum of the
  7 patches of each cluster, fires above a threshold;
* ``TriggerTDSCANEmulator``, the candidate firmware trigger: the traces are
  quantized to 4-bit codes and filtered by the TDSCAN kernel, in floating
  point or in fixed point (bit-exact with the HLS core), and it fires above a
  threshold.

Both take the R0 data of one telescope event and return the trigger output,
one boolean per patch and sample (432, T). They work on the 432 triplet
traces (the FADC output, one trace per patch of 3 pixels), which come from:

* real data: ``trigger_input_traces`` of the R0 data, computed by the FPGA itself, put in
  patch_sw_id order (see ``readout_to_patch_order``) and with the module swaps
  (``SST1MEventSource.swapped_modules``) undone as for the pixels;
* simulation: the raw waveform minus the simulated pedestal, through the FADC
  of the gateware (signed pixels, upper clip only, triplet sum clipped to 0-255).

Every camera-wide operation is a product with a sparse DigiCam matrix:
``patch_matrix`` (pixels -> patches), ``cluster_7_matrix`` (patch7 clusters)
and the neighbour matrix derived from it (TDSCAN).

The ``sst1mpipe-trigger`` tool (``sst1mpipe/scripts/sst1mpipe_trigger_tool.py``)
runs both emulators on a file. Configuration examples:
``sst1mpipe/resources/config/sst1mpipe_trigger_emulator_{mc,data}.json``.
"""
import numpy as np
from ctapipe.core import TelescopeComponent
from ctapipe.core.traits import CaselessStrEnum, Dict, Float, FloatTelescopeParameter, Int, List
from scipy.sparse import identity

from sst1mpipe.constants import PATCH_ID_INPUT_SORT_IDS, PATCH_ID_OUTPUT_SORT_IDS
from sst1mpipe.instrument.camera import DigiCam
from sst1mpipe.trigger import fixed_point

N_PIXELS = 1296
N_PATCHES = 432
PATCHES_PER_SECTOR = 144


# ---------------------------------------------------------------------------
# Camera trigger geometry
# ---------------------------------------------------------------------------

def neighbor_matrix():
    """(432, 432) sparse matrix, 1 where two patches are neighbours: the patch7 clusters without their centre."""
    return (DigiCam.cluster_7_matrix - identity(N_PATCHES, dtype=np.int64, format="csr")).tocsr()


# Order of the three camera sectors (micro_crate in camera_config.cfg) in the
# trigger readout, per telescope. Measured by comparing trigger_input_traces
# with the triplets rebuilt from the pixel waveforms (correlation 0.997 for
# tel 21, 0.999 for tel 22), and confirmed on trigger_output_patch7 (identical
# to the emulated patch7 output, sample by sample, on both telescopes). The
# sst1mpipe reader assumes (1, 2, 3).
READOUT_SECTOR_ORDER = {21: (1, 2, 3), 22: (2, 3, 1)}


def _camera_config_columns():
    """camera_config.cfg columns used for the readout order, per patch_sw_id."""
    with open(DigiCam.config_file) as f:
        rows = [line.split() for line in f if line.strip() and not line.startswith("#")]
    columns = {"micro_crate": 3, "patch_in_mod_fw": 13, "module_fw_id": 14, "patch_trig_out": 15}
    return {name: {int(row[11]): int(row[column]) for row in rows} for name, column in columns.items()}


def readout_to_patch_order(sector_order, signal="input"):
    """Indices putting a trigger signal, as read by sst1mpipe, in patch_sw_id order.

    ``signal="input"``: ``trigger_input_traces``; each sector fills a block of
    144 patches, ordered by module (module_fw_id), then by patch in the module
    (patch_in_mod_fw). ``signal="output"``: ``trigger_output_patch7``; each
    sector fills a block of 144 patches, ordered by patch_trig_out. Both from
    the official camera_config.cfg.
    """
    config = _camera_config_columns()
    sector = config["micro_crate"]
    if signal == "input":
        module_rank = {m: i for i, m in enumerate(sorted(set(config["module_fw_id"].values())))}
        in_sector = {p: 4 * module_rank[config["module_fw_id"][p]] + config["patch_in_mod_fw"][p] for p in sector}
        sort_ids = PATCH_ID_INPUT_SORT_IDS
    elif signal == "output":
        in_sector = config["patch_trig_out"]
        sort_ids = PATCH_ID_OUTPUT_SORT_IDS
    else:
        raise ValueError(f"signal must be 'input' or 'output', got {signal!r}")
    readout_position = np.array([
        PATCHES_PER_SECTOR * sector_order.index(sector[p]) + in_sector[p] for p in range(N_PATCHES)
    ])
    # The sst1mpipe reader puts readout position sort_ids[h] in row h.
    return np.argsort(sort_ids)[readout_position]


def module_swap_order(swapped_pixels):
    """Patch order that undoes the module swaps of ``SST1MEventSource.swapped_modules``.

    For some periods two modules were connected in each other's place. SST1MEventSource
    moves their pixels back, but the FPGA trigger traces stay as recorded:
    ``traces[order]`` puts each patch back on its pixels.
    ``swapped_pixels`` is the list of pixel id pairs of ``SST1MEventSource.swapped_pixels``.
    """
    # Same moves as SST1MEventSource: corrected pixel p shows recorded pixel pixel_order[p].
    pixel_order = np.arange(N_PIXELS)
    for pixels_1, pixels_2 in swapped_pixels:
        pixel_order[pixels_1], pixel_order[pixels_2] = pixel_order[pixels_2], pixel_order[pixels_1]

    # Row p of the CSR patch matrix holds the 3 pixels of patch p.
    triplets = DigiCam.patch_matrix.indices.reshape(N_PATCHES, 3)
    patch_of_pixel = np.empty(N_PIXELS, dtype=int)
    patch_of_pixel[triplets] = np.arange(N_PATCHES)[:, None]
    source_patches = patch_of_pixel[pixel_order[triplets]]  # (432, 3), one patch per row if whole modules move
    if np.any(source_patches != source_patches[:, :1]):
        raise ValueError("Swapped modules split a trigger patch, the trigger traces cannot follow the pixels")
    return source_patches[:, 0]


# ---------------------------------------------------------------------------
# Trigger stages
# ---------------------------------------------------------------------------

def fadc(waveform, baseline, patch_matrix=DigiCam.patch_matrix):
    """Triplet traces (432, T) from the raw waveforms (1296, T), like the FADC gateware.

    Each pixel: sample - floor(baseline), signed, clipped above at 4095 only.
    Each triplet: sum of its 3 pixels, clipped to [0, 255].
    """
    pixels = waveform.astype(np.int64) - np.floor(baseline).astype(np.int64)[:, None]
    pixels = np.minimum(pixels, 4095)
    return np.clip(patch_matrix.dot(pixels), 0, 255)


def score_quantizer(traces, edges):
    """Code of each sample = number of edges it reaches (value >= edge)."""
    return np.searchsorted(edges, traces, side="right")


# ---------------------------------------------------------------------------
# Emulators
# ---------------------------------------------------------------------------

class TriggerEmulator(TelescopeComponent):
    """Base class: R0 data of one telescope event -> trigger output (432, T).

    A child class defines ``threshold`` and ``scores(traces)``; the trigger
    fires on a patch and sample where the score is above the threshold of the
    telescope.
    """

    def __init__(self, subarray, config=None, parent=None, **kwargs):
        super().__init__(subarray=subarray, config=config, parent=parent, **kwargs)
        self._readout_to_patch = {tel: readout_to_patch_order(order) for tel, order in READOUT_SECTOR_ORDER.items()}
        self._module_swaps = {}  # tel_id -> patch order undoing the module swaps, set on the first event as sst1mpipe does

    def trigger_traces(self, event, tel_id, pedestal=None, swapped_pixels=()):
        """
        Triplet traces (432, T), the FADC output, in patch_sw_id order.

        Simulation: ``pedestal`` is the simulated pedestal per sample (n_channels, n_pixels), e.g.
        ``camera_monitorings[tel_id]["pedestal"] / camera_monitorings[tel_id]["n_ped_slices"]``
        of the simtel file. Real data: ``swapped_pixels`` are the pixel id pairs of the wrongly
        connected modules at the time of the event (``SST1MEventSource.swapped_pixels``).
        """
        if event.simulation is not None:
            waveform = event.r0.tel[tel_id].waveform[0]
            return fadc(waveform, pedestal[0])

        if tel_id not in self._readout_to_patch:
            raise ValueError(f"No trigger readout order known for telescope {tel_id}, see READOUT_SECTOR_ORDER")
        traces = np.asarray(event.r0.tel[tel_id].trigger_input_traces)
        if not np.all(np.isfinite(traces)):
            raise ValueError(
                f"Telescope {tel_id}: no trigger_input_traces in this file (the reader filled them with NaN), "
                "the trigger cannot be emulated on it"
            )
        return traces[self._readout_to_patch[tel_id]][self.module_swap(tel_id, swapped_pixels)].astype(np.int64)

    def module_swap(self, tel_id, swapped_pixels):
        """Patch order undoing the module swaps of this telescope (real data)."""
        if tel_id not in self._module_swaps:
            self._module_swaps[tel_id] = module_swap_order(swapped_pixels)
        return self._module_swaps[tel_id]

    def scores(self, traces):
        """Score of each patch and sample (432, T), compared to the threshold."""
        raise NotImplementedError

    def trigger_output(self, traces, tel_id):
        """Trigger output (432, T) bool and the scores it is decided on, from the triplet traces."""
        scores = self.scores(traces)
        return scores > self.threshold.tel[tel_id], scores

    def __call__(self, event, tel_id):
        """Trigger output (432, T) bool and the scores it is decided on, from the R0 data."""
        return self.trigger_output(self.trigger_traces(event, tel_id), tel_id)


class TriggerPatch7Emulator(TriggerEmulator):
    """patch7: sum of the 7 patches of each cluster (fewer at the camera edge)."""

    threshold = FloatTelescopeParameter(
        default_value=242.0, help="patch7 fires when a 7-patch cluster sum is above this value"
    ).tag(config=True)

    def scores(self, traces):
        return DigiCam.cluster_7_matrix.dot(traces)

    def hardware_output(self, event, tel_id, swapped_pixels=()):
        """patch7 output recorded by the camera (432, T) bool, in patch_sw_id order (real data)."""
        output = np.asarray(event.r0.tel[tel_id].trigger_output_patch7)
        if not np.all(np.isfinite(output)):
            raise ValueError(f"Telescope {tel_id}: no trigger_output_patch7 in this file")
        order = readout_to_patch_order(READOUT_SECTOR_ORDER[tel_id], signal="output")
        return output[order][self.module_swap(tel_id, swapped_pixels)].astype(bool)


class TriggerTDSCANEmulator(TriggerEmulator):
    """TDSCAN: 4-bit codes filtered in space (patch + neighbours) and time.

    For each patch and sample t, with L = 2 * eps_t + 1 time taps::

        conv[tap, t] = w_centre[tap] * code[patch, t'] + w_neighbours[tap] * sum of code[neighbour, t']
        score[t]     = sum over the taps of conv[tap, t],    t' = t + tap - eps_t

    Samples outside the readout window and neighbours outside the camera count
    as 0. ``quantize_step=None``: floating point, the weights are used as
    given. Otherwise fixed point: weights and codes are quantized, ``conv`` and
    ``score`` are requantized to the ``convolution_accumulator`` and
    ``temporal_accumulator`` formats, as in the firmware.
    """

    QUANTIZE_STEP_KEYS = (
        "input", "ring_weights", "convolution_accumulator", "convolution_rescale_shift",
        "temporal_accumulator", "temporal_rescale_shift",
    )

    threshold = FloatTelescopeParameter(
        default_value=9.03125, help="TDSCAN fires when a filtered score is above this value"
    ).tag(config=True)
    score_quantizer_edges = List(
        Float(), default_value=[16.0, 24.0, 32.0, 40.0, 48.0, 56.0, 64.0, 72.0, 80.0, 88.0, 96.0, 104.0, 112.0, 120.0, 128.0],
        help="Strictly increasing edges: the code of a sample is the number of edges it reaches",
    ).tag(config=True)
    eps_t = Int(2, help="Temporal half-width, the kernel spans 2 * eps_t + 1 samples").tag(config=True)
    ring_weights = List(
        List(Float()),
        default_value=[[0.5, -0.0078125], [-0.5, -0.03125], [-0.5, -0.125], [0.5, 0.0625], [0.5, 0.25]],
        help="Weights, one row per time tap, [centre, neighbours]",
    ).tag(config=True)
    quantize_step = Dict(
        default_value=None, allow_none=True,
        help="Fixed-point formats of the HLS core (see sst1mpipe_trigger_emulator_data.json), null for floating point",
    ).tag(config=True)
    overflow_mode = CaselessStrEnum(["AP_SAT", "AP_WRAP"], default_value="AP_SAT", help="Fixed point only").tag(config=True)
    quantization_mode = CaselessStrEnum(["AP_TRN", "AP_RND"], default_value="AP_TRN", help="Fixed point only").tag(config=True)

    def __init__(self, subarray, config=None, parent=None, **kwargs):
        super().__init__(subarray=subarray, config=config, parent=parent, **kwargs)
        self.neighbors = neighbor_matrix()

        self.edges = np.array(self.score_quantizer_edges, dtype=float)
        if np.any(np.diff(self.edges) <= 0):
            raise ValueError("TriggerTDSCANEmulator.score_quantizer_edges must be strictly increasing")

        weights = np.array(self.ring_weights, dtype=float)
        if weights.shape != (2 * self.eps_t + 1, 2):
            raise ValueError(
                f"ring_weights must have 2 * eps_t + 1 = {2 * self.eps_t + 1} rows of [centre, neighbours], "
                f"got shape {weights.shape}"
            )

        self.fixed_point = self.quantize_step is not None
        if not self.fixed_point:
            self.weights = weights
            return

        step = self.quantize_step
        missing = set(self.QUANTIZE_STEP_KEYS) - set(step)
        if missing:
            raise ValueError(f"quantize_step is missing {sorted(missing)}, expected all of {list(self.QUANTIZE_STEP_KEYS)}")
        for key in ("convolution_rescale_shift", "temporal_rescale_shift"):
            if step[key] < 0 or step[key] != int(step[key]):
                raise ValueError(f"quantize_step {key} must be an integer >= 0, got {step[key]}")

        self.input_format = fixed_point.parse(step["input"])
        weight_format = fixed_point.parse(step["ring_weights"])
        self.conv_format = fixed_point.parse(step["convolution_accumulator"])
        self.score_format = fixed_point.parse(step["temporal_accumulator"])
        self.conv_shift = int(step["convolution_rescale_shift"])
        self.score_shift = int(step["temporal_rescale_shift"])

        # Full-precision registers before each requantization: a sum of n terms
        # needs ceil(log2(n)) more integer bits (7 neighbours -> +3, 5 taps -> +3).
        self.conv_full_format = fixed_point.Format(
            signed=self.input_format.signed or weight_format.signed,
            int_bits=self.input_format.int_bits + weight_format.int_bits + int(np.ceil(np.log2(7))),
            frac_bits=self.input_format.frac_bits + weight_format.frac_bits,
        )
        self.score_full_format = fixed_point.Format(
            signed=self.conv_format.signed,
            int_bits=self.conv_format.int_bits + int(np.ceil(np.log2(len(weights)))),
            frac_bits=self.conv_format.frac_bits,
        )
        self.weights = fixed_point.to_codes(weights, weight_format, self.overflow_mode, self.quantization_mode)

    def scores(self, traces):
        return self.filter(score_quantizer(traces, self.edges))

    def filter(self, codes):
        """TDSCAN scores (432, T) in real units, for the score-quantizer codes (432, T)."""
        if self.fixed_point:
            codes = fixed_point.to_codes(codes, self.input_format, self.overflow_mode, self.quantization_mode)
        n_patches, n_samples = codes.shape

        # eps_t zero samples on each side; tap k looks at sample t + k - eps_t = index t + k.
        padded = np.zeros((n_patches, n_samples + 2 * self.eps_t), dtype=self.weights.dtype)
        padded[:, self.eps_t:self.eps_t + n_samples] = codes
        neighbor_sums = self.neighbors.dot(padded)

        score = np.zeros((n_patches, n_samples), dtype=self.weights.dtype)
        for tap, (w_centre, w_neighbors) in enumerate(self.weights):
            window = slice(tap, tap + n_samples)
            conv = w_centre * padded[:, window] + w_neighbors * neighbor_sums[:, window]
            if self.fixed_point:
                conv = fixed_point.requantize(conv, self.conv_full_format, self.conv_format, self.conv_shift, self.overflow_mode)
            score += conv

        if not self.fixed_point:
            return score
        score = fixed_point.requantize(score, self.score_full_format, self.score_format, self.score_shift, self.overflow_mode)
        return fixed_point.to_real(score, self.score_format)
