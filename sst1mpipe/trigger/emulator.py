"""Software emulation of the SST-1M DigiCam trigger.

For each telescope event, two triggers are computed from the 432 triplet
traces (the FADC output, one trace per patch of 3 pixels):

* patch7, the trigger running on the telescope: sum of the 7 patches of each
  cluster, fires when a cluster sum is above a threshold;
* TDSCAN, the candidate firmware trigger: the traces are quantized to 4-bit
  codes, filtered by the TDSCAN kernel, in floating point or in fixed point
  (bit-exact with the HLS core), and it fires when a filtered score is above
  a threshold.

Where the triplet traces come from:

* real data: ``trigger_input_traces``, computed by the FPGA itself, with the
  module swaps of ``inverted_module_list.json`` undone as for the pixels;
* simulation: the raw waveform minus the simulated pedestal, with the FADC of
  the gateware (signed pixels, upper clip only, triplet sum clipped to 0-255).

To enable it, copy the ``TriggerEmulator`` section of
``sst1mpipe/data/sst1mpipe_trigger_emulator_mc.json`` (simulations) or
``sst1mpipe_trigger_emulator_data.json`` (real data) into the config given to
``sst1mpipe_r0_dl1 --config``. Without this section the script is unchanged.
Thresholds are telescope parameters, e.g. ``[["type", "*", 350], ["id", 21, 225]]``.

Usage in ``sst1mpipe_r0_dl1``::

    emulator = TriggerEmulator(subarray=source.subarray, config=config)
    for event in source:
        emulator.process(event, event.trigger.tels_with_trigger, is_simulation, is_pedestal)  # before bad pixel removal
        ...
        if not emulator.keep(event, list(event.trigger.tels_with_trigger)):      # optional event filter
            continue
        ...
    emulator.write(output_file)
"""
import logging
from collections import Counter

import numpy as np
import tables
from astropy.table import Table
from ctapipe.containers import Map
from ctapipe.core import TelescopeComponent
from ctapipe.core.traits import Bool, CaselessStrEnum, Dict, Float, FloatTelescopeParameter, Int, List
from ctapipe.io import read_table

from sst1mpipe.constants import PATCH_ID_INPUT_SORT_IDS
from sst1mpipe.instrument.camera import DigiCam
from sst1mpipe.io.containers import R0CameraContainer, R0Container, SST1MArrayEventContainer, SST1MContainer
from sst1mpipe.trigger import fixed_point
from sst1mpipe.utils import get_swaped_modules, get_telescopes

N_PIXELS = 1296


# ---------------------------------------------------------------------------
# Camera trigger geometry
# ---------------------------------------------------------------------------

def read_trigger_geometry():
    """Patch tables of the trigger, from the DigiCam camera object (camera_config.cfg).

    Patches are in patch_sw_id order, pixels in pixel_sw_id order.

    Returns
    -------
    triplets: (432, 3) array
        The 3 pixels of each patch.
    clusters: list of 432 arrays
        The patches of the patch7 cluster centred on each patch (7, fewer at the edge).
    neighbors: (432, 7) array
        The same clusters as TDSCAN neighbourhoods: the patch itself in
        column 3, its neighbours in the other columns, -1 when missing at the
        camera edge.
    """
    triplets = np.array([np.nonzero(row)[0] for row in DigiCam.patch_matrix.toarray()])
    clusters = [np.nonzero(row)[0] for row in DigiCam.cluster_7_matrix.toarray()]

    neighbors = np.full((len(clusters), 7), -1)
    for patch, cluster in enumerate(clusters):
        others = [p for p in cluster if p != patch]
        neighbors[patch, 3] = patch
        neighbors[patch, [0, 1, 2, 4, 5, 6][:len(others)]] = others

    return triplets, clusters, neighbors


# Order of the three camera sectors (micro_crate in camera_config.cfg) in the
# readout of trigger_input_traces, per telescope. Measured by comparing
# trigger_input_traces with the triplets rebuilt from the pixel waveforms
# (correlation 0.997 for tel 21, 0.999 for tel 22). The sst1mpipe reader
# (PATCH_ID_INPUT) assumes (1, 2, 3).
READOUT_SECTOR_ORDER = {21: (1, 2, 3), 22: (2, 3, 1)}


def readout_to_patch_order(sector_order):
    """Indices putting ``trigger_input_traces`` (as read by sst1mpipe) in patch_sw_id order.

    In the readout, each sector fills a block of 144 patches, ordered by module
    (module_fw_id) then by patch in the module (patch_in_mod_fw), as given by
    the official camera_config.cfg.
    """
    # module_fw_id and patch_in_mod_fw are not kept by the DigiCam object: read its config file.
    with open(DigiCam.config_file) as f:
        rows = [line.split() for line in f if line.strip() and not line.startswith("#")]
    sector, module, patch_in_module = {}, {}, {}
    for row in rows:
        patch = int(row[11])                  # patch_sw_id
        sector[patch] = int(row[3])           # micro_crate
        module[patch] = int(row[14])          # module_fw_id
        patch_in_module[patch] = int(row[13]) # patch_in_mod_fw
    module_rank = {m: i for i, m in enumerate(sorted(set(module.values())))}
    readout_position = np.array([
        144 * sector_order.index(sector[p]) + 4 * module_rank[module[p]] + patch_in_module[p]
        for p in range(len(sector))
    ])
    # sst1mpipe row h holds readout position PATCH_ID_INPUT_SORT_IDS[h].
    row_of_readout_position = np.argsort(PATCH_ID_INPUT_SORT_IDS)
    return row_of_readout_position[readout_position]


# ---------------------------------------------------------------------------
# Trigger stages
# ---------------------------------------------------------------------------

def module_swap_order(triplets, swapped_modules):
    """Patch order that undoes the module swaps of ``inverted_module_list.json``.

    For some periods two modules were connected in each other's place. sst1mpipe
    moves their pixels back (``swap_modules_r0wf``), but the FPGA trigger traces
    stay as recorded: ``traces[order]`` puts each patch back on its pixels.
    ``swapped_modules`` is the list of pixel mask pairs of ``get_swaped_modules``.
    """
    # Same moves as swap_modules_r0wf: corrected pixel p shows recorded pixel pixel_order[p].
    pixel_order = np.arange(N_PIXELS)
    for mask_1, mask_2 in swapped_modules:
        pixels_1, pixels_2 = np.flatnonzero(mask_1), np.flatnonzero(mask_2)
        pixel_order[pixels_1], pixel_order[pixels_2] = pixel_order[pixels_2], pixel_order[pixels_1]

    patch_of_pixel = np.empty(N_PIXELS, dtype=int)
    patch_of_pixel[triplets.ravel()] = np.repeat(np.arange(len(triplets)), triplets.shape[1])
    source_patches = patch_of_pixel[pixel_order[triplets]]  # (n_patches, 3), one patch per row if whole modules move
    if np.any(source_patches != source_patches[:, :1]):
        raise ValueError("Swapped modules split a trigger patch, the trigger traces cannot follow the pixels")
    return source_patches[:, 0]


def as_sst1m_event(event):
    """Copy a simulated ctapipe array event into an SST1M one (shared fields, no data copy)."""
    sst1m_event = SST1MArrayEventContainer()
    for name in event.keys():
        setattr(sst1m_event, name, event[name])
    # Container defaults are shared between instances: give each event its own R0 telescope map,
    # so a telescope does not keep the trigger output of a previous event.
    sst1m_event.sst1m = SST1MContainer(r0=R0Container(tel=Map(R0CameraContainer)))
    return sst1m_event


def fadc(waveform, baseline, triplets):
    """Triplet traces (432, T) from the raw waveforms (1296, T), like the FADC gateware.

    Each pixel: sample - floor(baseline), signed, clipped above at 4095 only.
    Each triplet: sum of its 3 pixels, clipped to [0, 255].
    """
    pixels = waveform.astype(np.int64) - np.floor(baseline).astype(np.int64)[:, None]
    pixels = np.minimum(pixels, 4095)
    return np.clip(pixels[triplets].sum(axis=1), 0, 255)


def patch7(traces, clusters):
    """Cluster sums (432, T): for each patch, the sum of the 7 patches of its cluster."""
    return np.array([traces[cluster].sum(axis=0) for cluster in clusters])


def score_quantizer(traces, edges):
    """Code of each sample = number of edges it reaches (value >= edge)."""
    return np.searchsorted(edges, traces, side="right")


class TDSCAN:
    """TDSCAN filter, in floating point or in fixed point (bit-exact with the HLS core).

    For each patch and sample t, with L = 2 * eps_t + 1 time taps::

        conv[tap, patch, t] = sum over the 7 neighbours of weight[tap, ring] * code[neighbour, t + tap - eps_t]
        score[patch, t]     = sum over the taps of conv[tap, patch, t]

    Samples outside the readout window and neighbours outside the camera count as 0.

    ``quantize_step=None``: floating point, the weights are used as given.
    Otherwise fixed point: weights and codes are quantized, ``conv`` and
    ``score`` are requantized to the ``convolution_accumulator`` and
    ``temporal_accumulator`` formats, as in the firmware.
    """

    # Ring of each of the 7 neighbour slots of the neighbour table (0 = the patch itself, in slot 3).
    RINGS = np.array([1, 1, 1, 0, 1, 1, 1])
    QUANTIZE_STEP_KEYS = (
        "input", "ring_weights", "convolution_accumulator", "convolution_rescale_shift",
        "temporal_accumulator", "temporal_rescale_shift",
    )

    def __init__(self, neighbors, eps_t, ring_weights, quantize_step, overflow_mode, quantization_mode):
        if neighbors.shape[1] != 7 or not np.array_equal(neighbors[:, 3], np.arange(len(neighbors))):
            raise ValueError("The TDSCAN neighbour table must have 7 columns with the patch itself in column 3")
        self.neighbors = neighbors
        # Missing neighbours (-1, camera edge) point to an extra all-zero patch, added after the last one.
        self.neighbors_padded = np.where(neighbors < 0, len(neighbors), neighbors)
        self.eps_t = eps_t
        ring_weights = np.array(ring_weights, dtype=float)
        if ring_weights.shape != (2 * eps_t + 1, 2):
            raise ValueError(
                f"ring_weights must have 2 * eps_t + 1 = {2 * eps_t + 1} rows of [centre, neighbours], "
                f"got shape {ring_weights.shape}"
            )
        # Weight of each (tap, neighbour slot).
        weights = ring_weights[:, self.RINGS]

        self.fixed_point = quantize_step is not None
        if not self.fixed_point:
            self.weights = weights
            return

        missing = set(self.QUANTIZE_STEP_KEYS) - set(quantize_step)
        if missing:
            raise ValueError(f"quantize_step is missing {sorted(missing)}, expected all of {list(self.QUANTIZE_STEP_KEYS)}")

        self.overflow = overflow_mode
        self.quantization = quantization_mode
        self.input_format = fixed_point.parse(quantize_step["input"])
        weight_format = fixed_point.parse(quantize_step["ring_weights"])
        self.conv_format = fixed_point.parse(quantize_step["convolution_accumulator"])
        self.score_format = fixed_point.parse(quantize_step["temporal_accumulator"])
        for key in ("convolution_rescale_shift", "temporal_rescale_shift"):
            if quantize_step[key] < 0 or quantize_step[key] != int(quantize_step[key]):
                raise ValueError(f"quantize_step {key} must be an integer >= 0, got {quantize_step[key]}")
        self.conv_shift = int(quantize_step["convolution_rescale_shift"])
        self.score_shift = int(quantize_step["temporal_rescale_shift"])

        # Full-precision registers before each requantization: a sum of n terms
        # needs ceil(log2(n)) more integer bits (7 neighbours -> +3, 5 taps -> +3).
        n_taps = 2 * eps_t + 1
        self.conv_full_format = fixed_point.Format(
            signed=self.input_format.signed or weight_format.signed,
            int_bits=self.input_format.int_bits + weight_format.int_bits + int(np.ceil(np.log2(7))),
            frac_bits=self.input_format.frac_bits + weight_format.frac_bits,
        )
        self.score_full_format = fixed_point.Format(
            signed=self.conv_format.signed,
            int_bits=self.conv_format.int_bits + int(np.ceil(np.log2(n_taps))),
            frac_bits=self.conv_format.frac_bits,
        )

        self.weights = fixed_point.to_codes(weights, weight_format, overflow_mode, quantization_mode)

    def __call__(self, codes):
        """Scores (432, T) in real units, for the score-quantizer codes (432, T)."""
        if self.fixed_point:
            codes = fixed_point.to_codes(codes, self.input_format, self.overflow, self.quantization)
        n_patches, n_samples = codes.shape

        # Zero padding: eps_t samples on each side, and the extra all-zero patch.
        padded = np.zeros((n_patches + 1, n_samples + 2 * self.eps_t), dtype=self.weights.dtype)
        padded[:n_patches, self.eps_t:self.eps_t + n_samples] = codes
        neighbor_codes = padded[self.neighbors_padded]  # (432, 7, T + 2 eps_t)

        score = np.zeros((n_patches, n_samples), dtype=self.weights.dtype)
        for tap in range(len(self.weights)):
            # This tap looks at sample t + tap - eps_t, i.e. index t + tap of the padded trace.
            window = neighbor_codes[:, :, tap:tap + n_samples]
            conv = (window * self.weights[tap][None, :, None]).sum(axis=1)
            if self.fixed_point:
                conv = fixed_point.requantize(conv, self.conv_full_format, self.conv_format, self.conv_shift, self.overflow)
            score += conv

        if not self.fixed_point:
            return score
        score = fixed_point.requantize(score, self.score_full_format, self.score_format, self.score_shift, self.overflow)
        return fixed_point.to_real(score, self.score_format)


# ---------------------------------------------------------------------------
# Emulator used by sst1mpipe_r0_dl1
# ---------------------------------------------------------------------------

class TriggerEmulator(TelescopeComponent):
    """Runs patch7 and TDSCAN on every telescope event and records the results.

    Configured by the ``TriggerEmulator`` section of the sst1mpipe config. The
    default thresholds give 7 kHz (the camera readout limit) on the simulated
    medium NSB. On real data the NSB changes from night to night, so the
    thresholds have to be set for the analysed runs.
    """

    enabled = Bool(False, help="Emulate the trigger in sst1mpipe_r0_dl1").tag(config=True)
    filter_events = Bool(False, help="Write only the events where the filter_by trigger fired").tag(config=True)
    filter_by = CaselessStrEnum(["tdscan", "patch7"], default_value="tdscan", help="Trigger used by filter_events").tag(config=True)
    restrict_cleaning_to_tdscan_mask = Bool(
        False, help="Keep only the cleaned pixels where TDSCAN fired"
    ).tag(config=True)

    patch7_threshold = FloatTelescopeParameter(
        default_value=242.0, help="patch7 fires when a 7-patch cluster sum is above this value"
    ).tag(config=True)
    tdscan_threshold = FloatTelescopeParameter(
        default_value=9.03125, help="TDSCAN fires when a filtered score is above this value"
    ).tag(config=True)
    score_quantizer_edges = List(
        Float(), default_value=[16.0, 24.0, 32.0, 40.0, 48.0, 56.0, 64.0, 72.0, 80.0, 88.0, 96.0, 104.0, 112.0, 120.0, 128.0],
        help="Strictly increasing edges: the code of a sample is the number of edges it reaches",
    ).tag(config=True)
    eps_t = Int(2, help="TDSCAN temporal half-width, the kernel spans 2 * eps_t + 1 samples").tag(config=True)
    ring_weights = List(
        List(Float()),
        default_value=[[0.5, -0.0078125], [-0.5, -0.03125], [-0.5, -0.125], [0.5, 0.0625], [0.5, 0.25]],
        help="TDSCAN weights, one row per time tap, [centre, neighbours]",
    ).tag(config=True)
    quantize_step = Dict(
        default_value=None, allow_none=True,
        help="Fixed-point formats of the HLS core (see sst1mpipe_trigger_emulator_data.json), null for floating point",
    ).tag(config=True)
    overflow_mode = CaselessStrEnum(["AP_SAT", "AP_WRAP"], default_value="AP_SAT", help="Fixed point only").tag(config=True)
    quantization_mode = CaselessStrEnum(["AP_TRN", "AP_RND"], default_value="AP_TRN", help="Fixed point only").tag(config=True)

    def __init__(self, subarray, config=None, parent=None, **kwargs):
        super().__init__(subarray=subarray, config=config, parent=parent, **kwargs)
        self.triplets, self.clusters, neighbors = read_trigger_geometry()
        self.readout_to_patch = {tel: readout_to_patch_order(order) for tel, order in READOUT_SECTOR_ORDER.items()}
        self.module_swaps = {}  # tel_id -> patch order undoing the module swaps, set on the first event as sst1mpipe does

        self.score_edges = np.array(self.score_quantizer_edges, dtype=float)
        if np.any(np.diff(self.score_edges) <= 0):
            raise ValueError("TriggerEmulator.score_quantizer_edges must be strictly increasing")
        self.tdscan = TDSCAN(
            neighbors,
            eps_t=self.eps_t,
            ring_weights=self.ring_weights,
            quantize_step=self.quantize_step,
            overflow_mode=self.overflow_mode,
            quantization_mode=self.quantization_mode,
        )

        self.results = {}  # (obs_id, event_id, tel_id) -> results, for the events that can be written
        self.current = {}  # tel_id -> results of the event being processed
        self.current_is_pedestal = False
        self.counts = {}   # tel_id -> event counters, logged at the end

    def triplet_traces(self, event, tel_id, is_simulation):
        """FADC output (432, T), in patch_sw_id order."""
        if is_simulation:
            waveform = event.r0.tel[tel_id].waveform[0]
            baseline = event.mon.tel[tel_id].calibration.pedestal_per_sample[0]
            return fadc(waveform, baseline, self.triplets)
        if tel_id not in self.readout_to_patch:
            raise ValueError(f"No trigger readout order known for telescope {tel_id}, see READOUT_SECTOR_ORDER")
        traces = np.asarray(event.sst1m.r0.tel[tel_id].trigger_input_traces)
        if not np.all(np.isfinite(traces)):
            raise ValueError(
                f"Telescope {tel_id}: no trigger_input_traces in this file (the reader filled them with NaN), "
                "the trigger cannot be emulated on it"
            )
        if tel_id not in self.module_swaps:
            self.module_swaps[tel_id] = module_swap_order(self.triplets, get_swaped_modules(event))
        return traces[self.readout_to_patch[tel_id]][self.module_swaps[tel_id]]

    def run(self, traces, tel_id):
        """patch7 and TDSCAN on one set of triplet traces (432, T), with the thresholds of telescope ``tel_id``.

        Returns the summary and the TDSCAN binary output (432, T).
        """
        cluster_sums = patch7(traces, self.clusters)
        scores = self.tdscan(score_quantizer(traces, self.score_edges))
        tdscan_output = scores > self.tdscan_threshold.tel[tel_id]
        patch7_threshold = self.patch7_threshold.tel[tel_id]

        fired_patches = tdscan_output.any(axis=1)
        pixel_mask = np.zeros(N_PIXELS, dtype=bool)
        pixel_mask[self.triplets[fired_patches].ravel()] = True

        summary = {
            "patch7_max": int(cluster_sums.max()),
            "patch7": bool(cluster_sums.max() > patch7_threshold),
            "tdscan_max": float(scores.max()),
            "tdscan": bool(fired_patches.any()),
            "tdscan_pixel_mask": pixel_mask,
        }
        return summary, tdscan_output

    def process(self, event, tel_ids, is_simulation, is_pedestal):
        """Emulate the trigger of every telescope of the event. Call on raw (R0) data.

        Pedestal events (real data, random triggers) are emulated and counted,
        which tells how often the trigger fires on NSB alone, but they are not
        stored: they are never written in the DL1 file.

        The TDSCAN output is also set on ``event.sst1m.r0.tel[tel_id]``. The
        event must be an SST1M event (see ``as_sst1m_event``).
        """
        self.current = {}
        self.current_is_pedestal = is_pedestal
        for tel_id in tel_ids:
            result, tdscan_output = self.run(self.triplet_traces(event, tel_id, is_simulation), tel_id)
            self.current[tel_id] = result
            event.sst1m.r0.tel[tel_id].trigger_output_tdscan = tdscan_output

            kind = "pedestal" if is_pedestal else "shower"
            counts = self.counts.setdefault(tel_id, Counter())
            counts[kind] += 1
            counts[kind + "_tdscan"] += result["tdscan"]
            counts[kind + "_patch7"] += result["patch7"]

            if not is_pedestal:
                key = (event.index.obs_id, event.index.event_id, tel_id)
                if key in self.results:
                    raise ValueError(f"Event (obs_id, event_id, tel_id) = {key} seen twice: trigger results would be mixed up")
                self.results[key] = result

    def keep(self, event, tel_ids):
        """Apply the event filter. Returns False if no telescope fired.

        Telescopes that did not fire are removed from the event so they are
        neither used in the stereo reconstruction nor written.
        """
        if not self.filter_events:
            return True
        fired = [tel_id for tel_id in tel_ids if self.current[tel_id][self.filter_by]]
        for tel_id in tel_ids:
            if tel_id not in fired:
                event.dl1.tel.pop(tel_id, None)
                event.trigger.tel.pop(tel_id, None)
                del self.results[(event.index.obs_id, event.index.event_id, tel_id)]
        event.trigger.tels_with_trigger = np.array(fired, dtype=int)
        return len(fired) > 0

    def restrict_cleaning(self, image_processor):
        """Make the image cleaning keep only pixels where TDSCAN fired."""
        image_processor.clean = TriggerMaskedCleaner(image_processor.clean, self)

    def write(self, output_file):
        """Add the trigger results to the DL1 file.

        * columns ``trigger_emulation_*`` in ``/dl1/event/telescope/parameters/tel_XXX``;
        * the TDSCAN pixel masks in ``/dl1/event/telescope/trigger_emulation/tel_XXX``.
        """
        for tel in get_telescopes(output_file):
            tel_id = int(tel.split("_")[-1])
            path = "/dl1/event/telescope/parameters/" + tel
            params = read_table(output_file, path)
            keys = zip(params["obs_id"], params["event_id"], strict=True)
            rows = [self.results[(obs_id, event_id, tel_id)] for obs_id, event_id in keys]
            for name in ("tdscan_max", "tdscan", "patch7_max", "patch7"):
                params["trigger_emulation_" + name] = [row[name] for row in rows]
            params.write(output_file, path=path, overwrite=True, append=True)

            masks = Table({
                "obs_id": params["obs_id"],
                "event_id": params["event_id"],
                "tdscan_pixel_mask": np.array([row["tdscan_pixel_mask"] for row in rows]).reshape(-1, N_PIXELS),
            })
            masks.write(output_file, path="/dl1/event/telescope/trigger_emulation/" + tel, overwrite=True, append=True)

        with tables.open_file(output_file, mode="a") as f:
            f.root._v_attrs["trigger_emulation_filter_events"] = self.filter_events
            f.root._v_attrs["trigger_emulation_filter_by"] = self.filter_by

    def log_counts(self):
        for tel_id, counts in sorted(self.counts.items()):
            for kind in ("shower", "pedestal"):
                if counts[kind] > 0:
                    logging.info(
                        "Trigger emulation, tel %d, %s events: %d, TDSCAN fired %d, patch7 fired %d",
                        tel_id, kind, counts[kind], counts[kind + "_tdscan"], counts[kind + "_patch7"],
                    )
            if self.filter_events:
                n_kept = sum(1 for key in self.results if key[2] == tel_id)
                logging.info("Trigger emulation, tel %d: %d events kept by the %s filter", tel_id, n_kept, self.filter_by)


class TriggerMaskedCleaner:
    """Image cleaner wrapper: the cleaning mask is restricted to the TDSCAN pixel mask.

    Pedestal events are cleaned normally, so the pedestal monitoring of the
    r0->dl1 script is not affected. Every attribute is forwarded to the
    wrapped cleaner, so the script can keep setting
    ``image_processor.clean.average_charge`` etc.
    """

    def __init__(self, cleaner, emulator):
        object.__setattr__(self, "cleaner", cleaner)
        object.__setattr__(self, "emulator", emulator)

    def __call__(self, tel_id, image, arrival_times=None):
        mask = self.cleaner(tel_id, image, arrival_times)
        if self.emulator.current_is_pedestal:
            return mask
        mask = mask & self.emulator.current[tel_id]["tdscan_pixel_mask"]
        if mask.sum() <= 1:  # same guard as the sst1mpipe cleaners: timing fails with <= 1 pixel
            mask[:] = False
        return mask

    def __getattr__(self, name):
        return getattr(self.cleaner, name)

    def __setattr__(self, name, value):
        setattr(self.cleaner, name, value)
