import numpy as np
from astropy.time import Time

S_TO_NS = np.uint64(1e9)

def camera_clock_to_time(local_camera_clock):
    """
    Convert the camera clock (ns, TAI scale) to an astropy Time with ns precision, see
    https://github.com/cta-observatory/ctapipe_io_nectarcam/issues/24
    """
    localtime = np.uint64(local_camera_clock)
    full_seconds = localtime // S_TO_NS
    fractional_seconds = (localtime % S_TO_NS) / S_TO_NS
    return Time(full_seconds, fractional_seconds, format='unix_tai')


def local_time_to_time(local_time_second, local_time_nanosec):

    local_time_second = np.uint64(local_time_second)
    local_time_nanosec = np.uint64(local_time_nanosec)
    fractional_seconds = local_time_nanosec / S_TO_NS

    return Time(local_time_second, fractional_seconds, format='unix_tai')
