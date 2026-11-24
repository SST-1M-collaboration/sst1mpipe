"""
DL3 files of sst1mpipe read with gammapy: synthetic photon lists written by photon_df_to_fits,
indexed by create_hdu_index, and the theta2 analysis of get_theta2_from_dl3
"""
import astropy.units as u
import numpy as np
import pandas as pd
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from gammapy.data import DataStore
from gammapy.maps import MapAxis

from sst1mpipe.analysis.analysis import get_theta2_from_dl3
from sst1mpipe.io.dl2_to_dl3 import create_hdu_index, photon_df_to_fits
from sst1mpipe.io.sst1m_event_source import SUBARRAY_DESCRIPTION

TARGET = SkyCoord(83.6331, 22.0145, unit="deg")
# wobble: pointing 0.4 deg from the target
POINTING = TARGET.directional_offset_by(90 * u.deg, 0.4 * u.deg)
N_SIGNAL, N_BACKGROUND = 300, 3000
T_START, DURATION = 1768946400.0, 1200.0  # unix s
OBS_IDS = (1, 2)
# IRF of the photon lists (header IRF), looked up in $CALDB by create_hdu_index
IRF_DIR = "data/sst1m_stereo/test/bcf/rf_gc0.75"


def photon_list(seed=1):
    """N_SIGNAL photons around the target (0.03 deg), N_BACKGROUND photons uniform within 1.5 deg of the pointing"""
    rng = np.random.default_rng(seed)
    signal = TARGET.directional_offset_by(
        rng.uniform(0, 360, N_SIGNAL) * u.deg, np.abs(rng.normal(0, 0.03, N_SIGNAL)) * u.deg
    )
    background = POINTING.directional_offset_by(
        rng.uniform(0, 360, N_BACKGROUND) * u.deg, np.sqrt(rng.uniform(0, 1.5**2, N_BACKGROUND)) * u.deg
    )
    n = N_SIGNAL + N_BACKGROUND
    return pd.DataFrame({
        "tel_setup": "stereo",
        "local_time": np.sort(T_START + rng.uniform(0, DURATION, n)),
        "array_ra": POINTING.ra.deg,
        "array_dec": POINTING.dec.deg,
        "true_alt_tel": 70.0,
        "true_az_tel": 180.0,
        "reco_alt": 70.0,
        "reco_az": 180.0,
        "reco_ra": np.concatenate([signal.ra.deg, background.ra.deg]),
        "reco_dec": np.concatenate([signal.dec.deg, background.dec.deg]),
        "event_id": np.arange(n),
        "reco_energy": 1.0,
        "sst1mpipe_version": "test",
    })


@pytest.fixture(scope="module")
def dl3_dir(tmp_path_factory):
    """DL3 directory: two observations of the same photon list, with their index files"""
    directory = tmp_path_factory.mktemp("dl3")
    photons = photon_list()
    files = []
    for obs_id in OBS_IDS:
        hdus = photon_df_to_fits(
            photons, obs_id=obs_id, target_name="Crab", target_pos=TARGET,
            config={"analysis": {"global_gammaness_cut": 0.75}}, RF_used="rf",
            start_t=T_START, end_t=T_START + DURATION, subarray=SUBARRAY_DESCRIPTION,
        )
        files.append(str(directory / f"dl3_{obs_id}.fits"))
        fits.HDUList(hdus).writeto(files[-1])
    # placeholder IRF file: the theta2 analysis reads the events only
    caldb = directory / "caldb"
    (caldb / IRF_DIR).mkdir(parents=True)
    fits.PrimaryHDU().writeto(caldb / IRF_DIR / "irf.fits")
    create_hdu_index(files, out_dir=str(directory), irf_dir=str(caldb))
    return directory


def test_data_store(dl3_dir):
    data_store = DataStore.from_dir(dl3_dir)
    assert sorted(data_store.obs_table["OBS_ID"]) == list(OBS_IDS)
    np.testing.assert_allclose(data_store.obs_table["LIVETIME"], DURATION)
    hdu_classes = dict(zip(data_store.hdu_table["HDU_TYPE"], data_store.hdu_table["HDU_CLASS"], strict=False))
    assert hdu_classes["psf"] == "psf_table" and hdu_classes["bkg"] == "bkg_2d"
    observation = data_store.obs(OBS_IDS[0], required_irf="all-optional")
    assert len(observation.events.table) == N_SIGNAL + N_BACKGROUND
    assert observation.get_pointing_icrs(observation.tmid).separation(POINTING) < 1 * u.arcsec


def test_theta2_from_dl3(dl3_dir):
    theta2_axis = MapAxis.from_bounds(0, 0.5, nbin=10, interp="lin", unit="deg2")
    counts_on, counts_off, alpha, event_counts = get_theta2_from_dl3(
        dl3_dir, target_coords=TARGET, theta2_axis=theta2_axis, theta_cut=0.1 * u.deg,
    )
    n_obs = len(OBS_IDS)
    assert len(counts_on) == len(counts_off) == theta2_axis.nbin
    # first bin (theta < 0.22 deg): the signal and the background
    background_density = N_BACKGROUND / (np.pi * 1.5**2)  # per deg2
    expected_background = background_density * np.pi * 0.05
    assert counts_on[0] == pytest.approx(n_obs * (N_SIGNAL + expected_background), rel=0.1)
    # 5 OFF regions, background only in the first bin
    assert counts_off[0] / 5 == pytest.approx(n_obs * expected_background, rel=0.2)
    assert event_counts.N_excess == pytest.approx(n_obs * N_SIGNAL, rel=0.1)
    assert event_counts.significance_lima > 20
    assert 0.1 < alpha < 0.3
    assert event_counts.t_elapsed.to_value(u.s) == pytest.approx(n_obs * DURATION)

    # selection of the observations
    counts_on_1, *_ = get_theta2_from_dl3(
        dl3_dir, good_obsids=np.array(OBS_IDS[:1]), target_coords=TARGET,
        theta2_axis=theta2_axis, theta_cut=0.1 * u.deg,
    )
    np.testing.assert_array_equal(n_obs * np.asarray(counts_on_1), counts_on)
