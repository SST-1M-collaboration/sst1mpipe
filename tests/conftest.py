"""
Shared test data: the SST-1M runs of the test data and a reader of the events of the zfits files.
The constants and functions are imported by the test modules (``from conftest import ...``).
"""
import numpy as np
from astropy.utils import iers
from protozfits import File

from sst1mpipe.resources import TEST_DATA_DIR

# the tests do not download the IERS tables (Earth orientation): the alt/az of the pointing
# are computed with the tables of astropy, precise enough for the tests (< 1 arcmin), without
# depending on the network
iers.conf.auto_download = False

ZFITS_DIR = TEST_DATA_DIR / "zfits"
S_TO_NS = 1_000_000_000

# observation of Mrk 421 (wobble W1) by both telescopes, with Cherenkov and pedestal events:
# the run of tel 22 is taken during the run of tel 21. Number of events (and of pedestal events)
# used in the tests: the 150 first events of tel 21, all the events of tel 22
OBSERVATION_RUNS = {
    21: dict(path=ZFITS_DIR / "SST1M1_20260120_1179.fits.fz", n_events=150, n_pedestals=40),
    22: dict(path=ZFITS_DIR / "SST1M2_20260120_1102.fits.fz", n_events=35, n_pedestals=13),
}
# dark run of tel 22: pedestal events only
DARK_FILE = ZFITS_DIR / "SST1M2_20260119_0007.fits.fz"
# runs of the transition between two wobbles (night of 2026-01-21)
TRANSITION_FILES = {
    21: ZFITS_DIR / "SST1M1_20260121_0206.fits.fz",
    22: ZFITS_DIR / "SST1M2_20260121_0585.fits.fz",
}


def read_zfits_events(path, events=slice(None)):
    """
    SWAT event id, White Rabbit time (s, ns and total ns) and camera event type of the events of
    a zfits file, read with protozfits
    """
    with File(str(path)) as f:
        # the slice is limited to the events of the file
        rows = list(f.Events[slice(*events.indices(len(f.Events)))])
    time_s = np.array([row.local_time_sec for row in rows], dtype=np.int64)
    time_ns = np.array([row.local_time_nanosec for row in rows], dtype=np.int64)
    return dict(
        event_id=np.array([row.arrayEvtNum for row in rows], dtype=np.int64),
        time_s=time_s,
        time_ns=time_ns,
        time=time_s * S_TO_NS + time_ns,
        event_type=np.array([row.event_type for row in rows]),
    )


