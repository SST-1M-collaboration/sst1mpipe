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

* real data: ``trigger_input_traces``, computed by the FPGA itself;
* simulation: the raw waveform minus the simulated pedestal, with the FADC of
  the gateware (signed pixels, upper clip only, triplet sum clipped to 0-255).

To enable it, copy the ``TriggerEmulator`` section of
``sst1mpipe/data/sst1mpipe_trigger_emulator_mc.json`` (simulations) or
``sst1mpipe_trigger_emulator_data.json`` (real data) into the config given to
``sst1mpipe_r0_dl1 --config``. Without this section the script is unchanged.

Usage in ``sst1mpipe_r0_dl1``::

    emulator = TriggerEmulator(config["TriggerEmulator"])
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
from importlib.resources import files

import numpy as np
import tables
from astropy.table import Table
from ctapipe.io import read_table

from sst1mpipe.io.containers import SST1MArrayEventContainer
from sst1mpipe.trigger import fixed_point
from sst1mpipe.utils import get_telescopes

N_PIXELS = 1296


# ---------------------------------------------------------------------------
# Camera trigger geometry
# ---------------------------------------------------------------------------

def read_trigger_geometry():
    """Read the patch tables shipped in ``sst1mpipe/data``.

    All tables use the patch order of ``sst1m_trigger_patches.csv``.

    Returns
    -------
    triplets: (432, 3) array
        The 3 pixels of each patch.
    clusters: list of 432 arrays
        The 7 patches of the patch7 cluster centred on each patch.
    neighbors: (432, 7) array
        The hexagonal neighbours of each patch used by TDSCAN, -1 when the
        patch is at the camera edge.
    hardware_to_csv: (432,) array
        ``hardware_to_csv[h]`` is the patch of hardware patch ``h``.
    """
    data = files("sst1mpipe.data")

    # Each line: the 3 pixels of the patch, then the 3 pixels of each neighbour patch.
    rows = []
    with data.joinpath("sst1m_trigger_patches.csv").open() as f:
        for line in f:
            if line.strip():
                rows.append(np.array(line.split(","), dtype=int).reshape(-1, 3))
    triplets = np.array([row[0] for row in rows])
    patch_of_triplet = {tuple(triplet): patch for patch, triplet in enumerate(triplets)}
    clusters = [np.array([patch_of_triplet[tuple(t)] for t in row]) for row in rows]

    with data.joinpath("sst1m_trigger_neighbors_eps1.csv").open() as f:
        neighbors = np.loadtxt(f, delimiter=",", skiprows=1, dtype=int)
    with data.joinpath("sst1m_trigger_patch_hw_to_csv.txt").open() as f:
        hardware_to_csv = np.loadtxt(f, dtype=int)

    return triplets, clusters, neighbors, hardware_to_csv


# ---------------------------------------------------------------------------
# Trigger stages
# ---------------------------------------------------------------------------

def as_sst1m_event(event):
    """Copy a simulated ctapipe array event into an SST1M one (shared fields, no data copy)."""
    sst1m_event = SST1MArrayEventContainer()
    for name in event.keys():
        setattr(sst1m_event, name, event[name])
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

    def __init__(self, neighbors, eps_t, ring_weights, quantize_step, overflow_mode, quantization_mode):
        if neighbors.shape[1] != 7 or not np.array_equal(neighbors[:, 3], np.arange(len(neighbors))):
            raise ValueError("The TDSCAN neighbour table must have 7 columns with the patch itself in column 3")
        self.neighbors = neighbors
        # Missing neighbours (-1, camera edge) point to an extra all-zero patch, added after the last one.
        self.neighbors_padded = np.where(neighbors < 0, len(neighbors), neighbors)
        self.eps_t = eps_t
        # Weight of each (tap, neighbour slot).
        weights = np.array(ring_weights, dtype=float)[:, self.RINGS]

        self.fixed_point = quantize_step is not None
        if not self.fixed_point:
            self.weights = weights
            return

        self.overflow = overflow_mode
        self.quantization = quantization_mode
        self.input_format = fixed_point.parse(quantize_step["input"])
        weight_format = fixed_point.parse(quantize_step["ring_weights"])
        self.conv_format = fixed_point.parse(quantize_step["convolution_accumulator"])
        self.score_format = fixed_point.parse(quantize_step["temporal_accumulator"])
        self.conv_shift = quantize_step["convolution_rescale_shift"]
        self.score_shift = quantize_step["temporal_rescale_shift"]

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

class TriggerEmulator:
    """Runs patch7 and TDSCAN on every telescope event and records the results.

    ``config`` is the ``TriggerEmulator`` section of the sst1mpipe config.
    """

    def __init__(self, config):
        self.triplets, self.clusters, neighbors, hardware_to_csv = read_trigger_geometry()
        self.csv_to_hardware = np.argsort(hardware_to_csv)

        self.patch7_threshold = config["patch7"]["threshold"]
        tdscan = config["tdscan"]
        self.score_edges = np.array(tdscan["score_quantizer_edges"], dtype=float)
        if np.any(np.diff(self.score_edges) <= 0):
            raise ValueError("TriggerEmulator.tdscan.score_quantizer_edges must be strictly increasing")
        # No quantize_step (or null): floating-point TDSCAN, the overflow and
        # quantization modes are then not used.
        quantize_step = tdscan.get("quantize_step")
        self.tdscan = TDSCAN(
            neighbors,
            eps_t=tdscan["eps_t"],
            ring_weights=tdscan["ring_weights"],
            quantize_step=quantize_step,
            overflow_mode=tdscan["overflow_mode"] if quantize_step else None,
            quantization_mode=tdscan["quantization_mode"] if quantize_step else None,
        )
        self.tdscan_threshold = tdscan["threshold"]

        self.filter_events = config["filter_events"]
        self.filter_by = config["filter_by"]
        if self.filter_by not in ("tdscan", "patch7"):
            raise ValueError(f"TriggerEmulator.filter_by must be 'tdscan' or 'patch7', got {self.filter_by!r}")

        self.results = {}  # (obs_id, event_id, tel_id) -> results, for the events that can be written
        self.current = {}  # tel_id -> results of the event being processed
        self.current_is_pedestal = False
        self.counts = {}   # tel_id -> event counters, logged at the end

    def triplet_traces(self, event, tel_id, is_simulation):
        """FADC output (432, T), in the patch order of sst1m_trigger_patches.csv."""
        if is_simulation:
            waveform = event.r0.tel[tel_id].waveform[0]
            baseline = event.mon.tel[tel_id].calibration.pedestal_per_sample[0]
            return fadc(waveform, baseline, self.triplets)
        traces = np.asarray(event.sst1m.r0.tel[tel_id].trigger_input_traces)
        return traces[self.csv_to_hardware]

    def run(self, traces):
        """patch7 and TDSCAN on one set of triplet traces (432, T).

        Returns the summary and the TDSCAN binary output (432, T).
        """
        cluster_sums = patch7(traces, self.clusters)
        scores = self.tdscan(score_quantizer(traces, self.score_edges))
        tdscan_output = scores > self.tdscan_threshold

        fired_patches = tdscan_output.any(axis=1)
        pixel_mask = np.zeros(N_PIXELS, dtype=bool)
        pixel_mask[self.triplets[fired_patches].ravel()] = True

        summary = {
            "patch7_max": int(cluster_sums.max()),
            "patch7": bool(cluster_sums.max() > self.patch7_threshold),
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
            result, tdscan_output = self.run(self.triplet_traces(event, tel_id, is_simulation))
            self.current[tel_id] = result
            event.sst1m.r0.tel[tel_id].trigger_output_tdscan = tdscan_output

            kind = "pedestal" if is_pedestal else "shower"
            counts = self.counts.setdefault(tel_id, Counter())
            counts[kind] += 1
            counts[kind + "_tdscan"] += result["tdscan"]
            counts[kind + "_patch7"] += result["patch7"]

            if not is_pedestal:
                self.results[(event.index.obs_id, event.index.event_id, tel_id)] = result

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
