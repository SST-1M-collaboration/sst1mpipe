import json
from importlib.resources import files

import numpy as np
import pytest
from traitlets import TraitError
from traitlets.config import Config

from sst1mpipe.trigger import fixed_point
from sst1mpipe.trigger.emulator import (
    READOUT_SECTOR_ORDER,
    TriggerEmulator,
    fadc,
    module_swap_order,
    read_trigger_geometry,
    readout_to_patch_order,
    score_quantizer,
)
from sst1mpipe.utils import get_subarray


def emulator(quantize_step=None):
    """The TDSCAN wPow2 trigger, in fixed point (e.g. Q14, as the HLS core) or floating point (None)."""
    return TriggerEmulator(subarray=get_subarray(), config=Config({"TriggerEmulator": {"quantize_step": quantize_step}}))


Q14 = {
    "input": "UQ4.0",
    "ring_weights": "SQ1.7",
    "convolution_accumulator": "SQ5.2",
    "convolution_rescale_shift": 3,
    "temporal_accumulator": "SQ5.2",
    "temporal_rescale_shift": 3,
}


def test_example_configs_are_valid():
    for name in ("sst1mpipe_trigger_emulator_mc.json", "sst1mpipe_trigger_emulator_data.json"):
        with files("sst1mpipe.data").joinpath(name).open() as f:
            TriggerEmulator(subarray=get_subarray(), config=Config(json.load(f)))


def test_disabled_unless_enabled():
    assert not TriggerEmulator(subarray=get_subarray(), config=Config({})).enabled
    assert not TriggerEmulator(subarray=get_subarray(), config=Config({"TriggerEmulator": {"patch7_threshold": 300}})).enabled


def test_invalid_config_is_refused():
    invalid = (
        {"patch7": {"threshold": 350}},                        # old nested format
        {"eps_t": 1},                                          # 5 rows of ring_weights for 3 taps
        {"ring_weights": [[0.5, 0.1, 0.2]] * 5},               # rows must be [centre, neighbours]
        {"quantize_step": {"input": "UQ4.0"}},                 # incomplete fixed-point formats
        {"quantize_step": {**Q14, "convolution_rescale_shift": 3.5}},  # shifts are whole bit counts
        {"score_quantizer_edges": [16.0, 8.0]},                # edges not increasing
    )
    for section in invalid:
        with pytest.raises((ValueError, TraitError)):
            TriggerEmulator(subarray=get_subarray(), config=Config({"TriggerEmulator": section}))


def test_thresholds_per_telescope():
    with files("sst1mpipe.data").joinpath("sst1mpipe_trigger_emulator_data.json").open() as f:
        emu = TriggerEmulator(subarray=get_subarray(), config=Config(json.load(f)))
    assert emu.patch7_threshold.tel[21] == 225
    assert emu.patch7_threshold.tel[22] == 350


def test_fixed_point_formats():
    fmt = fixed_point.parse("SQ5.2")
    assert (fmt.bits, fmt.min_code, fmt.max_code) == (7, -64, 63)
    assert fixed_point.parse("UQ4.0").max_code == 15
    assert fixed_point.parse("sq5.2") == fmt

    # Truncation is a floor, also for negative values; saturation clips to the range.
    weights = fixed_point.to_codes([0.5, -0.0078, 5.0], fixed_point.parse("SQ1.7"), "AP_SAT", "AP_TRN")
    assert weights.tolist() == [64, -1, 127]

    # SQ8.7 -> SQ5.2 with shift 3: right shift by 15 - 7 - 3 = 5 bits, then saturate.
    source, target = fixed_point.parse("SQ8.7"), fixed_point.parse("SQ5.2")
    codes = fixed_point.requantize([1000, -1000, 5000], source, target, shift=3, overflow_mode="AP_SAT")
    assert codes.tolist() == [31, -32, 63]


def test_fadc_keeps_negative_pixels():
    triplets = np.array([[0, 1, 2]])
    baseline = np.array([10.7, 10.0, 10.0])
    waveform = np.array([[5, 10, 300], [13, 12, 300], [14, 1000, 300]])  # 3 pixels, 3 samples
    # sample 0: (5-10) + (13-10) + (14-10) = 2, the negative pixel is not clipped to 0
    # sample 1: (10-10) + (12-10) + (1000-10) = 992 -> 255
    # sample 2: (300-10) * 3 -> 255
    assert fadc(waveform, baseline, triplets).tolist() == [[2, 255, 255]]
    assert fadc(np.array([[0], [0], [0]]), baseline, triplets).tolist() == [[0]]


def test_score_quantizer():
    edges = np.array([16.0, 24.0])
    assert score_quantizer(np.array([15, 16, 23, 24, 100]), edges).tolist() == [0, 1, 1, 2, 2]


def test_trigger_geometry():
    triplets, clusters, neighbors = read_trigger_geometry()
    assert sorted(triplets.ravel().tolist()) == list(range(1296))
    for patch, cluster in enumerate(clusters):
        assert patch in cluster
        assert len(cluster) <= 7
        # TDSCAN neighbourhoods are the patch7 clusters, with the patch itself in column 3.
        assert set(neighbors[patch][neighbors[patch] >= 0]) == set(cluster)
    assert np.array_equal(neighbors[:, 3], np.arange(432))


def test_readout_order():
    # Sectors in the order sst1mpipe's reader assumes: nothing to reorder.
    assert np.array_equal(readout_to_patch_order((1, 2, 3)), np.arange(432))
    # Telescope 2: same layout inside each sector, sectors shifted by one.
    tel2 = readout_to_patch_order(READOUT_SECTOR_ORDER[22])
    assert sorted(tel2.tolist()) == list(range(432))
    assert not np.array_equal(tel2, np.arange(432))


def test_module_swap_moves_trigger_traces_with_pixels():
    from astropy.io import ascii

    from sst1mpipe.utils.utils import MAPPING_FILE_PATH

    # Telescope 2, 2023-09-01 to 2024-07-18: modules 59 and 88 connected in each other's place.
    pix_maps = ascii.read(MAPPING_FILE_PATH)
    masks = []
    for module in (59, 88):
        mask = np.zeros(1296, dtype=bool)
        mask[pix_maps[pix_maps["module"] == module]["pixel_sw_id"]] = True
        masks.append(mask)
    triplets = read_trigger_geometry()[0]
    order = module_swap_order(triplets, [masks])

    rng = np.random.default_rng(0)
    waveform, baseline = rng.integers(0, 100, (1296, 50)), np.full(1296, 50.0)
    corrected = waveform.copy()  # what swap_modules_r0wf does to the pixels
    corrected[masks[0]], corrected[masks[1]] = waveform[masks[1]], waveform[masks[0]]

    assert np.array_equal(fadc(waveform, baseline, triplets)[order], fadc(corrected, baseline, triplets))
    assert np.count_nonzero(order != np.arange(432)) == 8
    assert np.array_equal(module_swap_order(triplets, []), np.arange(432))


def test_tdscan_single_impulse():
    emu = emulator(Q14)
    codes = np.zeros((432, 50), dtype=int)
    patch, sample = 200, 20
    codes[patch, sample] = 1

    scores = emu.tdscan(codes)

    # Only one tap sees the impulse at each output sample: the score is the
    # centre weight of that tap (+-0.5, SQ1.7 code +-64) requantized to SQ5.2:
    # floor(+-64 / 32) = +-2 -> +-0.5. Tap k is seen at sample 20 + eps_t - k.
    assert scores[patch, 18:23].tolist() == [0.5, 0.5, -0.5, -0.5, 0.5]
    assert np.count_nonzero(scores[patch]) == 5

    # A neighbour sees the ring-1 weights: floor(code / 32) of -1, -4, -16, 8, 32.
    neighbor = emu.tdscan.neighbors[patch, 0]
    assert scores[neighbor, 18:23].tolist() == [0.25, 0.0, -0.25, -0.25, -0.25]


def test_tdscan_single_impulse_floating_point():
    emu = emulator(None)
    codes = np.zeros((432, 50), dtype=int)
    patch, sample = 200, 20
    codes[patch, sample] = 1

    scores = emu.tdscan(codes)

    # Without quantization the scores are the weights themselves, tap 4 first.
    assert scores[patch, 18:23].tolist() == [0.5, 0.5, -0.5, -0.5, 0.5]
    neighbor = emu.tdscan.neighbors[patch, 0]
    assert scores[neighbor, 18:23].tolist() == [0.25, 0.0625, -0.125, -0.03125, -0.0078125]


def test_as_sst1m_event_gives_each_event_its_own_r0():
    from ctapipe.containers import ArrayEventContainer

    from sst1mpipe.trigger import as_sst1m_event

    first, second = as_sst1m_event(ArrayEventContainer()), as_sst1m_event(ArrayEventContainer())
    first.sst1m.r0.tel[1].trigger_output_tdscan = np.ones((432, 50), dtype=bool)
    assert 1 not in second.sst1m.r0.tel


def test_missing_trigger_traces_are_refused():
    from sst1mpipe.io.containers import SST1MArrayEventContainer

    event = SST1MArrayEventContainer()
    event.sst1m.r0.tel[22].trigger_input_traces = np.zeros((432, 50)) * np.nan   # what the reader stores
    with pytest.raises(ValueError, match="no trigger_input_traces"):
        emulator().triplet_traces(event, 22, is_simulation=False)


def test_duplicate_event_is_refused():
    from sst1mpipe.io.containers import SST1MArrayEventContainer

    emu = emulator()
    traces = np.zeros((432, 50), dtype=np.int64)
    emu.triplet_traces = lambda event, tel_id, is_simulation: traces
    event = SST1MArrayEventContainer()
    event.index.obs_id, event.index.event_id = 1, 7
    emu.process(event, [22], is_simulation=False, is_pedestal=False)
    with pytest.raises(ValueError, match="seen twice"):
        emu.process(event, [22], is_simulation=False, is_pedestal=False)
