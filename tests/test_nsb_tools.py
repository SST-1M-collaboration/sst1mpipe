"""NSB level and voltage drop of the SiPMs (sst1mpipe.utils.nsb_tools, sst1mpipe-nsb)"""
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pytest
import tables
from astropy.table import QTable
from ctapipe.core import run_tool
from ctapipe.core.tool import ToolConfigurationError
from ctapipe.io import read_table
from traitlets.config import Config

from conftest import OBSERVATION_RUNS, WAVEFORM_CHUNK_ARGS, run_monitoring
from sst1mpipe.calib import R0R1Calibrator
from sst1mpipe.io import load_r0_pedestals
from sst1mpipe.scripts.sst1mpipe_nsb_tool import HDF5_TABLE, NSBTool
from sst1mpipe.utils import get_subarray
from sst1mpipe.utils.nsb_tools import (
    SiPMVoltageDropModel,
    mean_pedestal_std,
    pedestal_time,
    plot_nsb_rate,
    read_pedestal_tables,
)

# values of the former functions of NSB_tools (VAR_to_shift, VAR_to_Gdrop, VAR_to_Idrop, VAR_to_NSB,
# VAR_to_NSB_photon_rate) for these variances of the ADC samples of the pedestal events (ADC²)
REFERENCE = {
    21: dict(
        variance=[5.0, 20.0, 50.0, 70.0, 80.0],
        baseline_shift=[1.7503467977856118, 10.223413556597748, 35.84673919799134, 70.83016286701888, 185.3584096583282],
        gain_drop=[0.9914758110947841, 0.950211975979369, 0.8254263801057822, 0.6550571068376181, 0.09730454496394159],
        charge_drop=[0.9859272117458037, 0.9178037550049541, 0.7117922168481496, 0.4305254905491681, -0.49028161365295886],
        nsb_rate=[21.016611444500583, 128.08437203684935, 517.0017647842718, 1287.2407852776603, 22677.745780995534],
        nsb_photon_rate=[60.38539626209981, 378.87745946783684, 1712.9673725123184, 5595.926991665008, -12859.360081919142],
    ),
    22: dict(
        variance=[5.0, 50.0, 150.0, 290.0, 350.0],
        baseline_shift=[0.9276576674276822, 25.005349454276825, 83.8607133094047, 208.36439908093007, 280.3385847340131],
        gain_drop=[0.9989053639524353, 0.9704936876439534, 0.9010443582949025, 0.7541300090845025, 0.6692004700138645],
        charge_drop=[0.9982903269189308, 0.9539151409557678, 0.8454447053707671, 0.6159844124938458, 0.48333598833521385],
        nsb_rate=[9.879513060784914, 274.1020980181327, 990.1124620108923, 2939.3374043092467, 4456.550022221279],
        nsb_photon_rate=[28.244570666887515, 796.7595692902021, 3014.93147352771, 10281.53356042591, 17629.41273320626],
    ),
}


@pytest.fixture(scope="module")
def subarray():
    return get_subarray()


@pytest.fixture(scope="module")
def model(subarray):
    return SiPMVoltageDropModel(subarray=subarray)


@pytest.mark.parametrize("tel_id", sorted(REFERENCE))
def test_default_model(model, tel_id):
    """the default parameters give the values of the former NSB_tools functions"""
    reference = REFERENCE[tel_id]
    variance = np.array(reference["variance"])
    for name in ("baseline_shift", "gain_drop", "charge_drop"):
        np.testing.assert_allclose(getattr(model, name)(tel_id, variance), reference[name], rtol=1e-12)
    np.testing.assert_allclose(model.nsb_rate(tel_id, variance).to_value(u.MHz), reference["nsb_rate"], rtol=1e-12)
    np.testing.assert_allclose(model.nsb_photon_rate(tel_id, variance).to_value(u.MHz), reference["nsb_photon_rate"], rtol=1e-12)
    # a scalar variance (median of the camera) gives a scalar
    assert np.ndim(model.charge_drop(tel_id, variance[1])) == 0


def test_sample_duration(model):
    """the duration of a sample is given by the sampling rate of the cameras (250 MHz)"""
    for tel_id in (21, 22):
        assert u.isclose(model.sample_duration(tel_id), 4 * u.ns)


def test_configured_model(subarray):
    """all the parameters are configurable per telescope"""
    spline = dict(knots=[0, 0, 100, 100], coefficients=[0, 50, 0, 0], degree=1)
    config = Config({"SiPMVoltageDropModel": {
        "baseline_shift_spline": [["id", 21, spline], ["id", 22, spline]],
        "gain_drop_slope": [["id", 21, -0.01], ["id", 22, -0.02]],
        "pde_drop_slope": [["id", 21, -0.001], ["id", 22, 0.0]],
        "pulse_integral": [["id", 21, 20.0], ["id", 22, 25.0]],
        "photon_detection_efficiency": 0.5,
    }})
    model = SiPMVoltageDropModel(subarray=subarray, config=config)
    # linear spline: shift = variance / 2
    np.testing.assert_allclose(model.baseline_shift(21, [10, 40]), [5, 20])
    np.testing.assert_allclose(model.gain_drop(22, 40), 1 - 0.02 * 20)
    np.testing.assert_allclose(model.charge_drop(21, 40), 1 - 0.011 * 20)
    rate = 20 / (20.0 * (1 - 0.01 * 20)) / (4 * u.ns)
    assert u.isclose(model.nsb_rate(21, 40), rate)
    rate = 20 / (20.0 * (1 - 0.011 * 20)) / (4 * u.ns) / 0.5
    assert u.isclose(model.nsb_photon_rate(21, 40), rate)


def test_telescope_without_parameters(subarray):
    """the simulated telescopes (1, 2) have no parameters"""
    model = SiPMVoltageDropModel(subarray=get_subarray())
    with pytest.raises(KeyError):
        model.charge_drop(1, 10.0)
    config = Config({"SiPMVoltageDropModel": {"pde_drop_slope": [["type", "*", None], ["id", 21, -0.003]]}})
    model = SiPMVoltageDropModel(subarray=subarray, config=config)
    with pytest.raises(ValueError, match="pde_drop_slope for the telescope 22"):
        model.charge_drop(22, 10.0)


def test_voltage_drop_of_the_calibrator(subarray):
    """the voltage drop correction of the R0R1Calibrator uses the configured model"""
    pedestal_std = np.full((1, 1296), 5.0)
    pedestal_std[0, :10] = np.nan
    default = R0R1Calibrator(subarray=subarray)
    model = SiPMVoltageDropModel(subarray=subarray)
    np.testing.assert_allclose(default.voltage_drop(21, pedestal_std), model.charge_drop(21, 25.0))

    config = Config({
        "R0R1Calibrator": {"voltage_drop_correction": "pixelwise"},
        "SiPMVoltageDropModel": {"pde_drop_slope": [["id", 21, 0.0], ["id", 22, 0.0]]},
    })
    calibrator = R0R1Calibrator(subarray=subarray, config=config)
    np.testing.assert_allclose(calibrator.voltage_drop(21, pedestal_std)[0, 10:], model.gain_drop(21, 25.0))


def test_voltage_drop_model_configured_on_the_command_line(tmp_path):
    """the model of the calibrator of the tools is configured from the command line"""
    # telescope parameters on the command line: tuples
    slopes = '[("id", 21, -0.001), ("id", 22, -0.002)]'
    tool = run_monitoring(
        OBSERVATION_RUNS[22]["path"], tmp_path / "monitoring.h5", "--max-events=5",
        f"--SiPMVoltageDropModel.pde_drop_slope={slopes}",
    )
    assert tool.r0_r1_calibrator.voltage_drop_model.pde_drop_slope.tel[22] == -0.002
    assert "SiPMVoltageDropModel" in {cls.__name__ for cls in tool.classes}


@pytest.fixture(scope="module")
def monitoring_files(observation_monitoring, tmp_path_factory):
    """
    monitoring files of the observation of tel 22 and of its 20 first events, with the chunks of
    WAVEFORM_CHUNK_ARGS, and a file without pedestal table
    """
    directory = tmp_path_factory.mktemp("monitoring")
    first_events = directory / "first_events.monitoring.h5"
    run_monitoring(OBSERVATION_RUNS[22]["path"], first_events, "--max-events=20", *WAVEFORM_CHUNK_ARGS)
    empty = directory / "empty.h5"
    tables.open_file(empty, "w").close()
    return [observation_monitoring[22], first_events], empty


@pytest.mark.parametrize("last_chunk_only", [False, True])
def test_read_pedestal_tables(monitoring_files, last_chunk_only):
    paths, empty = monitoring_files
    pedestals = [load_r0_pedestals(path) for path in paths]
    assert len(pedestals[0]) > len(pedestals[1]) > 1

    # the tables of the files sorted by name; the files without pedestal table are skipped
    table = read_pedestal_tables([empty, *paths], last_chunk_only=last_chunk_only)
    expected = [p[-1:] if last_chunk_only else p for _, p in sorted(zip(map(str, paths), pedestals, strict=True))]
    assert len(table) == sum(len(p) for p in expected)
    np.testing.assert_array_equal(table["event_id_end"], np.concatenate([p["event_id_end"] for p in expected]))
    assert set(table["tel_id"]) == {22}
    assert np.all(np.isfinite(mean_pedestal_std(table)))

    assert read_pedestal_tables([empty], last_chunk_only=last_chunk_only) is None


def test_plot_nsb_rate(observation_monitoring, model):
    table = load_r0_pedestals(observation_monitoring[22])
    for photon_rate in (False, True):
        ax = plot_nsb_rate(table, 22, model=model, photon_rate=photon_rate, color="red")
        (line,) = ax.get_lines()
        variance = mean_pedestal_std(table) ** 2
        rate = model.nsb_photon_rate(22, variance) if photon_rate else model.nsb_rate(22, variance)
        np.testing.assert_allclose(line.get_ydata(), rate.to_value(u.MHz))
        assert "photon" in ax.get_ylabel() if photon_rate else "photon" not in ax.get_ylabel()
        plt.close(ax.figure)


def run_nsb_tool(*args):
    tool = NSBTool()
    run_tool(tool, argv=list(map(str, args)), raises=True)
    return tool


def test_nsb_tool(observation_monitoring, model, tmp_path):
    """NSB table of the chunks of pedestal events of the monitoring files of both telescopes, and plot"""
    output, plot = tmp_path / "nsb.ecsv", tmp_path / "nsb.png"
    run_nsb_tool("-i", observation_monitoring[21], "-i", observation_monitoring[22], "-o", output, "--plot", plot)
    assert plot.exists()

    table = QTable.read(output)
    assert set(table["tel_id"]) == {21, 22}
    for tel_id in (21, 22):
        pedestals = load_r0_pedestals(observation_monitoring[tel_id])
        rows = table[table["tel_id"] == tel_id]
        assert len(rows) == len(pedestals)
        variance = mean_pedestal_std(pedestals) ** 2
        np.testing.assert_allclose(rows["pedestal_variance"], variance)
        np.testing.assert_allclose(rows["charge_drop"], model.charge_drop(tel_id, variance))
        np.testing.assert_allclose(rows["nsb_rate"].to_value(u.MHz), model.nsb_rate(tel_id, variance).to_value(u.MHz))
        assert np.all(np.abs((rows["time"] - pedestal_time(pedestals)).to_value(u.s)) < 1e-6)
        assert np.all(rows["n_events"] == np.asarray(pedestals["n_events"]).reshape(len(pedestals), -1).max(axis=1))
    # NSB of the test runs: about 100 MHz per pixel
    assert np.all((20 * u.MHz < table["nsb_rate"]) & (table["nsb_rate"] < 1000 * u.MHz))


def test_nsb_tool_configuration(observation_monitoring, tmp_path):
    """the model is configured from the command line; HDF5 output"""
    outputs = {pde: tmp_path / f"nsb_{pde}.h5" for pde in (0.35, 0.7)}
    for pde, output in outputs.items():
        run_nsb_tool("-i", observation_monitoring[22], "-o", output, f"--SiPMVoltageDropModel.photon_detection_efficiency={pde}")
    default, configured = (read_table(outputs[pde], HDF5_TABLE) for pde in (0.35, 0.7))
    np.testing.assert_allclose(configured["nsb_photon_rate"], default["nsb_photon_rate"] / 2)
    np.testing.assert_allclose(configured["nsb_rate"], default["nsb_rate"])


def test_nsb_tool_errors(monitoring_files, tmp_path):
    paths, empty = monitoring_files
    with pytest.raises(ToolConfigurationError, match="No statistics"):
        run_nsb_tool("-i", empty, "-o", tmp_path / "nsb.ecsv")
    with pytest.raises(ToolConfigurationError, match="format"):
        run_nsb_tool("-i", paths[0], "-o", tmp_path / "nsb.txt")
    existing = tmp_path / "existing.ecsv"
    existing.touch()
    with pytest.raises(ToolConfigurationError, match="exists"):
        run_nsb_tool("-i", paths[0], "-o", existing)
    run_nsb_tool("-i", paths[0], "-o", existing, "--overwrite")
    assert len(QTable.read(existing)) > 0
