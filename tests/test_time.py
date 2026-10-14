import numpy as np

from sst1mpipe.time import camera_clock_to_time, local_time_to_time

def test_time():

    local_time_nanosec = 383583536
    local_time_sec = 1769015262

    local_camera_clock = np.int64(local_time_sec * 1E9) + local_time_nanosec

    time_1 = camera_clock_to_time(local_camera_clock)
    time_2 = local_time_to_time(local_time_sec, local_time_nanosec)
    time_1.precision = 9

    assert time_1.iso == '2026-01-21 17:07:42.383583536'
    assert time_1 == time_2
