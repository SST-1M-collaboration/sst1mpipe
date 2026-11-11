import json

import numpy as np
import pytest
import tables
from ctapipe.core import run_tool
from scipy.sparse import csr_matrix
from traitlets import TraitError
from traitlets.config import Config

from conftest import OBSERVATION_RUNS, TRANSITION_FILES
from sst1mpipe.instrument.camera import DigiCam
from sst1mpipe.resources import CONFIG_DIR
from sst1mpipe.trigger import TriggerPatch7Emulator, TriggerTDSCANEmulator, fixed_point
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.trigger.emulator import BASELINE_SAMPLE, TriggerEmulator, score_quantizer
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


def test_compute_trigger_traces_keeps_negative_pixels():
    patch_matrix = csr_matrix(np.ones((1, 3), dtype=np.int64))  # one patch of 3 pixels
    baseline = np.array([10.7, 10.0, 10.0])
    waveform = np.array([[5, 10, 300], [13, 12, 300], [14, 1000, 300]])  # 3 pixels, 3 samples
    # sample 0: (5-10) + (13-10) + (14-10) = 2, the negative pixel is not clipped to 0
    # sample 1: (10-10) + (12-10) + (1000-10) = 992 -> 255
    # sample 2: (300-10) * 3 -> 255
    assert TriggerEmulator._compute_trigger_traces(waveform, baseline, patch_matrix=patch_matrix).tolist() == [[2, 255, 255]]
    assert TriggerEmulator._compute_trigger_traces(np.zeros((3, 1)), baseline, patch_matrix=patch_matrix).tolist() == [[0]]
    # pixels without data (pixel_flags 0) do not contribute: (5-10) + (13-10) = -2 -> 0, 0 + 2 = 2, 290 + 290 -> 255
    assert TriggerEmulator._compute_trigger_traces(waveform, baseline, pixel_flags=[1, 1, 0], patch_matrix=patch_matrix).tolist() == [[0, 2, 255]]


def test_score_quantizer():
    edges = np.array([16.0, 24.0])
    assert score_quantizer(np.array([15, 16, 23, 24, 100]), edges).tolist() == [0, 1, 1, 2, 2]


def test_trigger_geometry():
    # Every pixel in exactly one patch of 3, each CSR row listing its 3 pixels.
    assert np.array_equal(DigiCam.patch_matrix.indptr, np.arange(0, 1297, 3))
    assert sorted(DigiCam.patch_matrix.indices.tolist()) == list(range(1296))

    neighbors = DigiCam.neighbor_patch_matrix
    assert (neighbors != neighbors.T).nnz == 0                  # symmetric
    assert neighbors.diagonal().sum() == 0                      # a patch is not its own neighbour
    assert set(np.asarray(neighbors.sum(axis=1)).ravel()) <= {3, 4, 6}  # 6 inside, fewer at the edge
    assert (neighbors + np.eye(432, dtype=np.int64) != DigiCam.cluster_7_matrix).sum() == 0


@pytest.mark.parametrize("tel_id", sorted(TRANSITION_FILES))
def test_trigger_input_traces_match_the_waveforms(tel_id):
    """
    the trigger input traces of the camera (FPGA) are the patches rebuilt from the pixel waveforms:
    checks the sector order of the telescope and the patch order of the reader. The FPGA updates its
    baseline continuously: the traces are identical at the sample of the stored baseline, and close
    elsewhere (mean squared error 4.4 ADC^2 per trace for tel 21, 6.3 for tel 22, ~5000 to 10000
    with a wrong order of the traces)
    """
    squared_error, n_traces = 0, 0
    for event in SST1MEventSource(input_url=TRANSITION_FILES[tel_id]):
        r0 = event.r0.tel[tel_id]
        camera = np.asarray(r0.trigger_input_traces).astype(np.int64)
        rebuilt = TriggerEmulator._compute_trigger_traces(r0.waveform[0], r0.pedestal, pixel_flags=r0.pixel_flags)
        np.testing.assert_array_equal(camera[:, BASELINE_SAMPLE], rebuilt[:, BASELINE_SAMPLE])
        squared_error += np.sum((camera - rebuilt) ** 2)
        n_traces += len(camera)
    assert squared_error / n_traces < 10


# patch7 threshold of the camera in the test files: the only threshold for which the emulated patch7
# output of the trigger_input_traces is the camera output, on all the events of the file
PATCH7_THRESHOLD_OF_FILES = [
    (21, TRANSITION_FILES[21], 210),
    (22, TRANSITION_FILES[22], 280),
    (21, OBSERVATION_RUNS[21]["path"], 185),
    (22, OBSERVATION_RUNS[22]["path"], 240),
]


@pytest.mark.parametrize("tel_id, path, threshold", PATCH7_THRESHOLD_OF_FILES, ids=lambda v: getattr(v, "name", v))
def test_patch7_output_from_camera_input(tel_id, path, threshold):
    """the patch7 output emulated from the trigger input traces of the camera is the camera output"""
    emu = TriggerPatch7Emulator(subarray=get_subarray(), threshold=threshold)
    n_events = 0
    for event in SST1MEventSource(input_url=path, max_events=500):
        r0 = event.r0.tel[tel_id]
        output, _ = emu.trigger_output(np.asarray(r0.trigger_input_traces).astype(np.int64), tel_id)
        np.testing.assert_array_equal(output, np.asarray(r0.trigger_output_patch7).astype(bool))
        n_events += 1
    assert n_events > 0


def test_patch7_sums_clusters():
    emu = TriggerPatch7Emulator(subarray=get_subarray(), threshold=242)
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
    neighbor = DigiCam.neighbor_patch_matrix[patch].indices[0]
    assert scores[neighbor, 18:23].tolist() == [0.25, 0.0, -0.25, -0.25, -0.25]


def test_tdscan_single_impulse_floating_point():
    emu = tdscan()
    codes = np.zeros((432, 50), dtype=int)
    patch, sample = 200, 20
    codes[patch, sample] = 1

    scores = emu.filter(codes)

    # Without quantization the scores are the weights themselves, tap 4 first.
    assert scores[patch, 18:23].tolist() == [0.5, 0.5, -0.5, -0.5, 0.5]
    neighbor = DigiCam.neighbor_patch_matrix[patch].indices[0]
    assert scores[neighbor, 18:23].tolist() == [0.25, 0.0625, -0.125, -0.03125, -0.0078125]


def test_missing_trigger_traces_are_refused():
    from sst1mpipe.io.containers import SST1MArrayEventContainer

    event = SST1MArrayEventContainer()
    event.r0.tel[22].trigger_input_traces = None   # what the reader stores
    with pytest.raises(ValueError, match="no trigger_input_traces"):
        tdscan().trigger_traces(event, 22)


@pytest.mark.parametrize("traces", [None, np.empty(0, dtype=np.uint8)])
def test_trigger_traces_not_in_the_camera_event(traces):
    """The trigger traces missing in a camera event (None or empty field) are None in the R0 data"""
    from sst1mpipe.io.sst1m_event_source import prepare_trigger_input, read_trigger_traces

    assert read_trigger_traces(traces, prepare_trigger_input, "trigger_input_traces", 50, np.arange(432)) is None


def test_trigger_tool(tmp_path):
    from sst1mpipe.scripts.sst1mpipe_trigger_tool import TriggerTool

    config = tmp_path / "config.json"
    config.write_text(json.dumps(example_config("data")))
    output = tmp_path / "trigger.h5"
    input_file = TRANSITION_FILES[21]
    argv = ["-i", str(input_file), "-o", str(output), "-m", "20", "--config", str(config),
            "--save-outputs", "--compare-camera"]
    assert run_tool(TriggerTool(), argv=argv, cwd=tmp_path) == 0

    with tables.open_file(output) as f:
        table = f.root.trigger_emulation.tel_021
        assert len(table) > 0
        assert table.coldtypes["trigger_traces"] == np.dtype((np.uint8, (432, 50)))
        assert table.coldtypes["tdscan_output"] == np.dtype((bool, (432, 50)))
        assert {"patch7", "tdscan", "camera_patch7", "camera_patch7_mismatches"} <= set(table.colnames)
