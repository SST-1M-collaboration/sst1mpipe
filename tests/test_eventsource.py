from importlib.resources import files

from sst1mpipe.io.sst1m_event_source import SST1MEventSource

FILE_TEL_1 = str(files('sst1mpipe.resources.zfits').joinpath('SST1M1_20260121_0001.fits.fz'))
FILE_TEL_2 = str(files('sst1mpipe.resources.zfits').joinpath('SST1M2_20260121_0001.fits.fz'))

MAX_ITERATIONS = 10

def test_read_events():

    source = SST1MEventSource([FILE_TEL_1], max_events=MAX_ITERATIONS)

    i = 0
    for _ in source:
        i += 1
    assert i == MAX_ITERATIONS
