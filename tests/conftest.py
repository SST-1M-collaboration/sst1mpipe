from astropy.utils import iers

# the tests do not download the IERS tables (Earth orientation): the alt/az of the pointing
# are computed with the tables of astropy, precise enough for the tests (< 1 arcmin), without
# depending on the network
iers.conf.auto_download = False
