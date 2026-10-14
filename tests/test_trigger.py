import json
from importlib.resources import files

import numpy as np
import pytest
import tables
from ctapipe.core import run_tool
from scipy.sparse import csr_matrix
from traitlets import TraitError
from traitlets.config import Config

from sst1mpipe.instrument.camera import DigiCam
from sst1mpipe.resources import CONFIG_DIR
from sst1mpipe.trigger import TriggerPatch7Emulator, TriggerTDSCANEmulator, fixed_point
from sst1mpipe.trigger.emulator import (
    READOUT_SECTOR_ORDER,
    fadc,
    module_swap_order,
    neighbor_matrix,
    readout_to_patch_order,
    score_quantizer,
)
from sst1mpipe.utils import get_subarray

Q14 = {
    "input": "UQ4.0",
    "ring_weights": "SQ1.7",
    "convolution_accumulator": "SQ5.2",
    "convolution_rescale_shift": 3,
    "temporal_accumulator": "SQ5.2",
    "temporal_rescale_shift": 3,
}


def example_config(name):
    with (CONFIG_DIR / f"sst1mpipe_trigger_emulator_{name}.json").open() as f:
        return Config(json.load(f))


def tdscan(**section):
    """TDSCAN wPow2 emulator, floating point unless quantize_step is given."""
    return TriggerTDSCANEmulator(subarray=get_subarray(), config=Config({"TriggerTDSCANEmulator": section}))


def test_example_configs_are_valid():
    for name in ("mc", "data"):
        TriggerPatch7Emulator(subarray=get_subarray(), config=example_config(name))
        TriggerTDSCANEmulator(subarray=get_subarray(), config=example_config(name))


def test_invalid_tdscan_config_is_refused():
    invalid = (
        {"eps_t": 1},                                          # 5 rows of ring_weights for 3 taps
        {"ring_weights": [[0.5, 0.1, 0.2]] * 5},               # rows must be [centre, neighbours]
        {"quantize_step": {"input": "UQ4.0"}},                 # incomplete fixed-point formats
        {"quantize_step": {**Q14, "convolution_rescale_shift": 3.5}},  # shifts are whole bit counts
        {"score_quantizer_edges": [16.0, 8.0]},                # edges not increasing
    )
    for section in invalid:
        with pytest.raises((ValueError, TraitError)):
            tdscan(**section)


def test_thresholds_per_telescope():
    emu = TriggerPatch7Emulator(subarray=get_subarray(), config=example_config("data"))
    assert emu.threshold.tel[21] == 225
    assert emu.threshold.tel[22] == 350


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
    patch_matrix = csr_matrix(np.ones((1, 3), dtype=np.int64))  # one patch of 3 pixels
    baseline = np.array([10.7, 10.0, 10.0])
    waveform = np.array([[5, 10, 300], [13, 12, 300], [14, 1000, 300]])  # 3 pixels, 3 samples
    # sample 0: (5-10) + (13-10) + (14-10) = 2, the negative pixel is not clipped to 0
    # sample 1: (10-10) + (12-10) + (1000-10) = 992 -> 255
    # sample 2: (300-10) * 3 -> 255
    assert fadc(waveform, baseline, patch_matrix).tolist() == [[2, 255, 255]]
    assert fadc(np.zeros((3, 1)), baseline, patch_matrix).tolist() == [[0]]


def test_score_quantizer():
    edges = np.array([16.0, 24.0])
    assert score_quantizer(np.array([15, 16, 23, 24, 100]), edges).tolist() == [0, 1, 1, 2, 2]


def test_trigger_geometry():
    # Every pixel in exactly one patch of 3, each CSR row listing its 3 pixels.
    assert np.array_equal(DigiCam.patch_matrix.indptr, np.arange(0, 1297, 3))
    assert sorted(DigiCam.patch_matrix.indices.tolist()) == list(range(1296))

    neighbors = neighbor_matrix()
    assert (neighbors != neighbors.T).nnz == 0                  # symmetric
    assert neighbors.diagonal().sum() == 0                      # a patch is not its own neighbour
    assert set(np.asarray(neighbors.sum(axis=1)).ravel()) <= {3, 4, 6}  # 6 inside, fewer at the edge
    assert (neighbors + np.eye(432, dtype=np.int64) != DigiCam.cluster_7_matrix).sum() == 0


def test_readout_order():
    for signal in ("input", "output"):
        # Sectors in the order sst1mpipe's reader assumes: nothing to reorder.
        assert np.array_equal(readout_to_patch_order((1, 2, 3), signal), np.arange(432))
        # Telescope 2: same layout inside each sector, sectors shifted by one.
        tel2 = readout_to_patch_order(READOUT_SECTOR_ORDER[22], signal)
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
    order = module_swap_order([masks])

    rng = np.random.default_rng(0)
    waveform, baseline = rng.integers(0, 100, (1296, 50)), np.full(1296, 50.0)
    corrected = waveform.copy()  # what swap_modules_r0wf does to the pixels
    corrected[masks[0]], corrected[masks[1]] = waveform[masks[1]], waveform[masks[0]]

    assert np.array_equal(fadc(waveform, baseline)[order], fadc(corrected, baseline))
    assert np.count_nonzero(order != np.arange(432)) == 8
    assert np.array_equal(module_swap_order([]), np.arange(432))


def test_patch7_sums_clusters():
    emu = TriggerPatch7Emulator(subarray=get_subarray())
    traces = np.zeros((432, 50), dtype=np.int64)
    patch, sample = 200, 20
    traces[patch, sample] = 100

    output, scores = emu.trigger_output(traces, tel_id=21)

    # Every cluster containing the patch sees it: the patch and its 6 neighbours.
    cluster = DigiCam.cluster_7_matrix[patch].indices
    assert sorted(np.flatnonzero(scores[:, sample]).tolist()) == sorted(cluster.tolist())
    assert scores.sum() == 100 * len(cluster)
    assert not output.any()                                   # 100 < 242
    traces[patch, sample] = 250
    assert emu.trigger_output(traces, tel_id=21)[0][:, sample].sum() == len(cluster)


def test_tdscan_single_impulse():
    emu = tdscan(quantize_step=Q14)
    codes = np.zeros((432, 50), dtype=int)
    patch, sample = 200, 20
    codes[patch, sample] = 1

    scores = emu.filter(codes)

    # Only one tap sees the impulse at each output sample: the score is the
    # centre weight of that tap (+-0.5, SQ1.7 code +-64) requantized to SQ5.2:
    # floor(+-64 / 32) = +-2 -> +-0.5. Tap k is seen at sample 20 + eps_t - k.
    assert scores[patch, 18:23].tolist() == [0.5, 0.5, -0.5, -0.5, 0.5]
    assert np.count_nonzero(scores[patch]) == 5

    # A neighbour sees the ring-1 weights: floor(code / 32) of -1, -4, -16, 8, 32.
    neighbor = neighbor_matrix()[patch].indices[0]
    assert scores[neighbor, 18:23].tolist() == [0.25, 0.0, -0.25, -0.25, -0.25]


def test_tdscan_single_impulse_floating_point():
    emu = tdscan()
    codes = np.zeros((432, 50), dtype=int)
    patch, sample = 200, 20
    codes[patch, sample] = 1

    scores = emu.filter(codes)

    # Without quantization the scores are the weights themselves, tap 4 first.
    assert scores[patch, 18:23].tolist() == [0.5, 0.5, -0.5, -0.5, 0.5]
    neighbor = neighbor_matrix()[patch].indices[0]
    assert scores[neighbor, 18:23].tolist() == [0.25, 0.0625, -0.125, -0.03125, -0.0078125]


def test_missing_trigger_traces_are_refused():
    from sst1mpipe.io.containers import SST1MArrayEventContainer

    event = SST1MArrayEventContainer()
    event.sst1m.r0.tel[22].trigger_input_traces = np.zeros((432, 50)) * np.nan   # what the reader stores
    with pytest.raises(ValueError, match="no trigger_input_traces"):
        tdscan().trigger_traces(event, 22)


def test_trigger_tool(tmp_path):
    from sst1mpipe.scripts.sst1mpipe_trigger_tool import TriggerTool

    config = tmp_path / "config.json"
    config.write_text(json.dumps(example_config("data")))
    output = tmp_path / "trigger.h5"
    input_file = files("sst1mpipe.resources.zfits").joinpath("SST1M1_20260121_0001.fits.fz")
    argv = ["-i", str(input_file), "-o", str(output), "-m", "20", "--config", str(config),
            "--save-outputs", "--compare-camera"]
    assert run_tool(TriggerTool(), argv=argv, cwd=tmp_path) == 0

    with tables.open_file(output) as f:
        table = f.root.trigger_emulation.tel_021
        assert len(table) > 0
        assert table.coldtypes["trigger_traces"] == np.dtype((np.uint8, (432, 50)))
        assert table.coldtypes["tdscan_output"] == np.dtype((bool, (432, 50)))
        assert {"patch7", "tdscan", "camera_patch7", "camera_patch7_mismatches"} <= set(table.colnames)
