import warnings
from collections import Counter

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import angular_separation
from ctapipe.containers import EventType
from ctapipe.core import run_tool
from ctapipe.io import EventSource, SimTelEventSource, read_table

from sst1mpipe.io import get_dl1_info
from sst1mpipe.resources import MC_CONFIG_FILES, TEST_DATA_DIR
from sst1mpipe.scripts.sst1mpipe_process_tool import ProcessorTool

SIMTEL_DIR = TEST_DATA_DIR / "simtel"
GAMMA_ID, PROTON_ID = 0, 101
TEL_IDS = [1, 2]
N_PIXELS = 1296
N_SAMPLES = 50
POINTING_ALT = 70 * u.deg
POINTING_AZ = 175.196 * u.deg

# content of the simulation files (sim_telarray, NSBmed4, zenith 20 deg)
FILES = {
    "gamma_diffuse": dict(
        name="gamma_200_800E3GeV_20_20deg_ATM52_102168.corsika.gz.NSBmed4.simtel.gz",
        run=102168, primary=GAMMA_ID, n_events=376, energy_range=(0.2, 800) * u.TeV,
        viewcone=10 * u.deg, n_showers=3000, scatter_range=928 * u.m,
    ),
    "gamma_point": dict(
        name="gamma_point_200_800E3GeV_20_20deg_ATM52_300669.corsika.gz.NSBmed4.simtel.gz",
        run=300669, primary=GAMMA_ID, n_events=563, energy_range=(0.2, 800) * u.TeV,
        viewcone=0 * u.deg, n_showers=1000, scatter_range=928 * u.m,
    ),
    "proton": dict(
        name="proton_400_1300E3GeV_20_20deg_ATM52_204508.corsika.gz.NSBmed4.simtel.gz",
        run=204508, primary=PROTON_ID, n_events=626, energy_range=(0.4, 1300) * u.TeV,
        viewcone=10 * u.deg, n_showers=2000, scatter_range=1032 * u.m,
    ),
}


@pytest.fixture(scope="module", params=sorted(FILES))
def simtel_file(request):
    return FILES[request.param] | {"path": SIMTEL_DIR / FILES[request.param]["name"]}


@pytest.fixture(scope="module")
def events(simtel_file):
    """summary of all the events of the file"""
    summary = dict(
        obs_ids=set(), event_ids=[], primaries=Counter(), energies=[], directions=[],
        multiplicity=Counter(), pointing=set(), r0=set(), r1=set(), true_images=[],
        event_types=Counter(),
    )
    with SimTelEventSource(simtel_file["path"]) as source:
        for event in source:
            shower = event.simulation.shower
            summary["obs_ids"].add(event.index.obs_id)
            summary["event_ids"].append(event.index.event_id)
            summary["event_types"][event.trigger.event_type] += 1
            summary["primaries"][int(shower.shower_primary_id)] += 1
            summary["energies"].append(shower.energy.to_value(u.TeV))
            summary["directions"].append((shower.alt.to_value(u.deg), shower.az.to_value(u.deg)))
            summary["multiplicity"][len(event.trigger.tels_with_trigger)] += 1
            for tel_id in event.trigger.tels_with_trigger:
                pointing = event.pointing.tel[tel_id]
                summary["pointing"].add((round(pointing.altitude.to_value(u.deg), 3), round(pointing.azimuth.to_value(u.deg), 3)))
                summary["r0"].add((event.r0.tel[tel_id].waveform.shape, event.r0.tel[tel_id].waveform.dtype))
                summary["r1"].add((event.r1.tel[tel_id].waveform.shape, event.r1.tel[tel_id].waveform.dtype))
                summary["true_images"].append(event.simulation.tel[tel_id].true_image)
    return summary


def test_simtel_file_is_read_by_simtel_event_source(simtel_file):
    assert SimTelEventSource.is_compatible(simtel_file["path"])
    with EventSource(simtel_file["path"]) as source:
        assert isinstance(source, SimTelEventSource)
        assert source.is_simulation
        assert source.obs_ids == [simtel_file["run"]]


def test_subarray(simtel_file):
    with SimTelEventSource(simtel_file["path"]) as source:
        subarray = source.subarray

    assert sorted(subarray.tel_ids) == TEL_IDS
    for tel_id in TEL_IDS:
        tel = subarray.tel[tel_id]
        assert tel.camera.name == "DigiCam"
        assert tel.camera.geometry.n_pixels == N_PIXELS
        assert tel.camera.readout.n_channels == 1
        assert tel.camera.readout.n_samples == N_SAMPLES
        assert tel.camera.readout.sampling_rate.to_value(u.GHz) == pytest.approx(0.25)
        assert tel.optics.equivalent_focal_length.to_value(u.m) == pytest.approx(5.6)
    # the two telescopes are ~155 m apart
    distance = np.linalg.norm(subarray.positions[1] - subarray.positions[2])
    assert distance.to_value(u.m) == pytest.approx(155.3, abs=0.1)


def test_simulation_config(simtel_file):
    with SimTelEventSource(simtel_file["path"]) as source:
        config = source.simulation_config[simtel_file["run"]]
        assert source.atmosphere_density_profile is not None

    assert config.run_number == simtel_file["run"]
    assert config.atmosphere == 52
    e_min, e_max = simtel_file["energy_range"]
    assert config.energy_range_min.to_value(u.TeV) == pytest.approx(e_min.to_value(u.TeV))
    assert config.energy_range_max.to_value(u.TeV) == pytest.approx(e_max.to_value(u.TeV))
    assert config.spectral_index == -2
    assert config.max_viewcone_radius.to_value(u.deg) == pytest.approx(simtel_file["viewcone"].to_value(u.deg))
    assert config.min_viewcone_radius.to_value(u.deg) == 0
    assert config.n_showers == simtel_file["n_showers"]
    assert config.shower_reuse == 20
    assert config.max_scatter_range.to_value(u.m) == pytest.approx(simtel_file["scatter_range"].to_value(u.m))


def test_events(simtel_file, events):
    n_events = simtel_file["n_events"]

    assert events["obs_ids"] == {simtel_file["run"]}
    assert len(events["event_ids"]) == n_events
    assert len(set(events["event_ids"])) == n_events
    assert events["event_types"] == {EventType.SUBARRAY: n_events}
    assert events["primaries"] == {simtel_file["primary"]: n_events}
    # mono and stereo events
    assert set(events["multiplicity"]) == {1, 2}

    energies = np.array(events["energies"]) * u.TeV
    e_min, e_max = simtel_file["energy_range"]
    assert np.all((energies >= e_min) & (energies <= e_max))


def test_pointing_and_shower_directions(simtel_file, events):
    # both telescopes point to the same direction, zenith 20 deg
    assert events["pointing"] == {(POINTING_ALT.to_value(u.deg), POINTING_AZ.to_value(u.deg))}

    alt, az = np.array(events["directions"]).T * u.deg
    offset = angular_separation(az, alt, POINTING_AZ, POINTING_ALT)
    viewcone = simtel_file["viewcone"]
    if viewcone == 0:
        # point source: all the showers come from the pointing direction
        np.testing.assert_allclose(offset.to_value(u.deg), 0, atol=1e-3)
    else:
        # diffuse: within the view cone, not all from the center
        assert np.all(offset <= viewcone + 0.01 * u.deg)
        assert np.max(offset) > 1 * u.deg


def test_waveforms_and_true_images(events):
    assert events["r0"] == {((1, N_PIXELS, N_SAMPLES), np.dtype("uint16"))}
    assert events["r1"] == {((1, N_PIXELS, N_SAMPLES), np.dtype("float32"))}

    true_images = np.array(events["true_images"])
    assert true_images.shape[1] == N_PIXELS
    assert np.issubdtype(true_images.dtype, np.integer)
    assert np.all(true_images >= 0)
    # Cherenkov photo-electrons in most of the triggered images
    assert np.mean(true_images.sum(axis=1) > 0) > 0.9


def test_process_simtel_file(simtel_file, tmp_path):
    """sst1mpipe-process from R0 (simulated waveforms) up to the DL1 parameters"""
    output = tmp_path / "events.dl1.h5"
    run_tool(ProcessorTool(), argv=[
        f"--input={simtel_file['path']}",
        f"--output={output}",
        f"--config={MC_CONFIG_FILES['low']}",
    ], raises=True)

    n_events = simtel_file["n_events"]
    trigger = read_table(output, "/dl1/event/subarray/trigger")
    shower = read_table(output, "/simulation/event/subarray/shower")
    assert len(trigger) == n_events
    assert len(shower) == n_events
    assert set(shower["true_shower_primary_id"]) == {simtel_file["primary"]}

    tel_trigger = read_table(output, "/dl1/event/telescope/trigger")
    info = get_dl1_info(output)
    for tel_id in TEL_IDS:
        tel = f"tel_{tel_id:03d}"
        parameters = read_table(output, f"/dl1/event/telescope/parameters/{tel}")
        n_triggered = np.sum(tel_trigger["tel_id"] == tel_id)
        assert len(parameters) == n_triggered
        assert info[f"n_triggered_tel{tel_id}"][0] == n_triggered
        # showers survive the cleaning
        intensity = parameters["camera_frame_hillas_intensity"]
        n_images = np.isfinite(intensity).sum()
        assert n_images > 0
        assert np.all(intensity[np.isfinite(intensity)] > 0)
        # the fraction of images surviving the cleaning depends on the configuration: only reported
        if n_images < 0.2 * n_triggered:
            warnings.warn(
                f"{simtel_file['name']}, tel {tel_id}: only {n_images} of the {n_triggered} images"
                " survive the cleaning", stacklevel=1,
            )
        # true images and parameters of the simulation
        assert len(read_table(output, f"/simulation/event/telescope/images/{tel}")) == n_triggered
        assert len(read_table(output, f"/simulation/event/telescope/parameters/{tel}")) == n_triggered
        # impact distance of the stereo reconstruction (ShowerProcessor, DL2)
        impact = read_table(output, f"/dl2/event/telescope/impact/HillasReconstructor/{tel}")
        assert len(impact) == n_triggered

    # shower geometry reconstructed for the stereo events
    geometry = read_table(output, "/dl2/event/subarray/geometry/HillasReconstructor")
    assert len(geometry) == n_events
    multiplicity = np.bincount(tel_trigger["event_id"].astype(np.int64))[geometry["event_id"]]
    is_valid = geometry["HillasReconstructor_is_valid"]
    assert is_valid.sum() > 0
    assert np.all(multiplicity[is_valid] == 2)

    # distribution of all the simulated showers, triggered or not
    distribution = read_table(output, "/simulation/service/shower_distribution")
    assert len(distribution) == 1
    assert distribution["n_entries"][0] > n_events
