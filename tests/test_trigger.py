import json
from importlib.resources import files

import numpy as np

from sst1mpipe.trigger import fixed_point
from sst1mpipe.trigger.emulator import TriggerEmulator, fadc, read_trigger_geometry, score_quantizer


def trigger_config(quantize_step):
    """The TDSCAN wPow2 trigger, in fixed point (Q14, as the HLS core) or floating point (None)."""
    return {
        "enabled": True,
        "filter_events": False,
        "filter_by": "tdscan",
        "restrict_cleaning_to_tdscan_mask": False,
        "patch7": {"threshold": 222},
        "tdscan": {
            "score_quantizer_edges": [16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128],
            "eps_t": 2,
            "ring_weights": [[0.5, -0.0078125], [-0.5, -0.03125], [-0.5, -0.125], [0.5, 0.0625], [0.5, 0.25]],
            "quantize_step": quantize_step,
            "overflow_mode": "AP_SAT",
            "quantization_mode": "AP_TRN",
            "threshold": 8.109374,
        },
    }


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
            TriggerEmulator(json.load(f)["TriggerEmulator"])


def test_fixed_point_formats():
    fmt = fixed_point.parse("SQ5.2")
    assert (fmt.bits, fmt.min_code, fmt.max_code) == (7, -64, 63)
    assert fixed_point.parse("UQ4.0").max_code == 15

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
    triplets, clusters, neighbors, hardware_to_csv = read_trigger_geometry()
    assert sorted(triplets.ravel().tolist()) == list(range(1296))
    for patch, cluster in enumerate(clusters):
        assert cluster[0] == patch
        assert len(cluster) <= 7
    assert np.array_equal(neighbors[:, 3], np.arange(432))
    assert sorted(hardware_to_csv.tolist()) == list(range(432))


def test_tdscan_single_impulse():
    emulator = TriggerEmulator(trigger_config(Q14))
    codes = np.zeros((432, 50), dtype=int)
    patch, sample = 200, 20
    codes[patch, sample] = 1

    scores = emulator.tdscan(codes)

    # Only one tap sees the impulse at each output sample: the score is the
    # centre weight of that tap (+-0.5, SQ1.7 code +-64) requantized to SQ5.2:
    # floor(+-64 / 32) = +-2 -> +-0.5. Tap k is seen at sample 20 + eps_t - k.
    assert scores[patch, 18:23].tolist() == [0.5, 0.5, -0.5, -0.5, 0.5]
    assert np.count_nonzero(scores[patch]) == 5

    # A neighbour sees the ring-1 weights: floor(code / 32) of -1, -4, -16, 8, 32.
    neighbor = emulator.tdscan.neighbors[patch, 0]
    assert scores[neighbor, 18:23].tolist() == [0.25, 0.0, -0.25, -0.25, -0.25]


def test_tdscan_single_impulse_floating_point():
    emulator = TriggerEmulator(trigger_config(None))
    codes = np.zeros((432, 50), dtype=int)
    patch, sample = 200, 20
    codes[patch, sample] = 1

    scores = emulator.tdscan(codes)

    # Without quantization the scores are the weights themselves, tap 4 first.
    assert scores[patch, 18:23].tolist() == [0.5, 0.5, -0.5, -0.5, 0.5]
    neighbor = emulator.tdscan.neighbors[patch, 0]
    assert scores[neighbor, 18:23].tolist() == [0.25, 0.0625, -0.125, -0.03125, -0.0078125]
