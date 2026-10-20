"""
SST1MMonitoringSource (monitoring files of sst1mpipe-monitoring)
"""
import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import AltAz
from ctapipe.io import MonitoringSource, read_table
from ctapipe.io.hdf5dataformat import DL1_PIXEL_STATISTICS_GROUP
from ctapipe.io.monitoringtypes import MonitoringType

from conftest import OBSERVATION_RUNS
from sst1mpipe.io import R0_PIXEL_STATISTICS_GROUP
from sst1mpipe.io.monitoring_source import SST1MMonitoringSource, statistics_container
from sst1mpipe.io.sst1m_event_source import SST1MEventSource

# observation of Mrk 421 by tel 22
TEL_ID = 22
INPUT_FILE = OBSERVATION_RUNS[TEL_ID]["path"]
N_EVENTS, N_PEDESTALS = OBSERVATION_RUNS[TEL_ID]["n_events"], OBSERVATION_RUNS[TEL_ID]["n_pedestals"]


@pytest.fixture(scope="module")
def monitoring(observation_monitoring):
    """monitoring file of the observation of tel 22, rolling chunks of 5 pedestal events"""
    return observation_monitoring[TEL_ID]


def pointing_altaz(source, time):
    """alt/az of the pointing of the run (RA/Dec of the event source) at the site of the array"""
    return source.pointing.transform_to(AltAz(obstime=time, location=source.subarray.reference_location))


def test_monitoring_source_plugin(monitoring):
    with SST1MEventSource(INPUT_FILE) as source:
        monitoring_source = MonitoringSource.from_name(
            "SST1MMonitoringSource", subarray=source.subarray, input_files=[monitoring],
        )
        # the telescope of the run
        assert source.subarray.tel_ids.tolist() == monitoring_source.subarray.tel_ids.tolist() == [TEL_ID]
        assert isinstance(monitoring_source, SST1MMonitoringSource)


def test_fill_monitoring_container(monitoring):
    waveform = read_table(monitoring, f"{R0_PIXEL_STATISTICS_GROUP}/sky_pedestal_waveform/tel_{TEL_ID:03d}")
    dead_pixels = np.all(waveform["outlier_mask"][:, 0], axis=0)
    assert dead_pixels.sum() > 0
    with SST1MEventSource(INPUT_FILE) as source:
        monitoring_source = SST1MMonitoringSource(subarray=source.subarray, input_files=[monitoring])
        for event in source:
            # the event source gives the RA/Dec of the pointing only
            assert np.isnan(event.monitoring.tel[TEL_ID].pointing.altitude)
            altaz = pointing_altaz(source, event.trigger.time)
            monitoring_source.fill_monitoring_container(event)
            monitoring = event.monitoring.tel[TEL_ID]
            statistics = monitoring.camera.pixel_statistics

            # statistics of the ADC samples, (n_channels, n_pixels): the outliers of all the chunks
            # (dead pixels) are NaN
            std = statistics.pedestal_waveform.std
            assert std.shape == (1, 1296)
            assert np.all(np.isnan(std[0, dead_pixels]))
            np.testing.assert_allclose(np.nanmedian(std), np.nanmedian(waveform["std"]), rtol=0.05)
            # statistics of the images, channel 0
            assert statistics.pedestal_image.std.shape == (1296, )
            # pointing interpolated from the pointing table
            assert u.isclose(monitoring.pointing.altitude, altaz.alt, atol=0.01 * u.deg)
            assert u.isclose(monitoring.pointing.azimuth, altaz.az, atol=0.01 * u.deg)


def test_statistics_of_the_last_pedestal_events(monitoring):
    """
    the statistics at the time of an event are those of the last chunk of pedestal events ended
    before it (rolling chunks: the last pedestal events), the first chunk before its end
    """
    waveform = read_table(monitoring, f"{R0_PIXEL_STATISTICS_GROUP}/sky_pedestal_waveform/tel_{TEL_ID:03d}")
    assert len(waveform) > 2
    std = np.where(waveform["outlier_mask"], np.nan, waveform["std"])
    chunks_used = set()
    with SST1MEventSource(INPUT_FILE) as source:
        monitoring_source = SST1MMonitoringSource(subarray=source.subarray, input_files=[monitoring])
        for event in source:
            monitoring_source.fill_monitoring_container(event)
            ended = np.flatnonzero(waveform["time_end"] <= event.trigger.time)
            chunk = ended[-1] if len(ended) else 0
            chunks_used.add(chunk)
            np.testing.assert_array_equal(event.monitoring.tel[TEL_ID].camera.pixel_statistics.pedestal_waveform.std, std[chunk])
    # all the chunks are used
    assert chunks_used == set(range(len(waveform)))


def test_no_statistics_out_of_the_chunks(monitoring):
    """no statistics without chunk of pedestal events at the time of the event (timestamp_tolerance)"""
    monitoring_source = SST1MMonitoringSource(input_files=[monitoring], timestamp_tolerance=0 * u.s)
    image = read_table(monitoring, f"{DL1_PIXEL_STATISTICS_GROUP}/sky_pedestal_image/tel_{TEL_ID:03d}")

    def statistics(time):
        values = monitoring_source.get_values(
            MonitoringType.PIXEL_STATISTICS, time=time, tel_id=TEL_ID, subtype="sky_pedestal_image",
        )
        return statistics_container(values, channel=0)

    assert statistics(image["time_start"][0]).std.shape == (1296, )
    # the time scale of the time does not matter (TAI in the file)
    assert statistics(image["time_start"][0].utc).std.shape == (1296, )
    assert statistics(image["time_start"][0] - 1 * u.s).std is None
    assert statistics(image["time_end"][-1] + 1 * u.s).std is None


