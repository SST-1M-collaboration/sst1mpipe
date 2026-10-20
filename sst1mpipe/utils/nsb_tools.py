"""
NSB level of the SST-1M cameras and voltage drop of their SiPMs, from the pedestal events.

The NSB light shifts the baseline of the SiPMs and increases the variance of their ADC samples.
The NSB current causes a voltage drop on the bias resistor of the SiPMs, which reduces their gain
and their PDE linearly with the baseline shift (fig. 11.1 of "SiPM behaviour under continuous
light", https://arxiv.org/abs/1910.00348).

`SiPMVoltageDropModel` gives, from the variance of the ADC samples of the pedestal events (ADC²,
the square of their std, see sst1mpipe-monitoring):

- the baseline shift (ADC), with a spline fitted for each telescope;
- the gain drop and the charge drop (gain and PDE), the factors by which the gain and the charges
  are reduced (the voltage drop correction of the `~sst1mpipe.calib.R0R1Calibrator` divides the
  charges by the charge drop);
- the NSB rate per pixel, in p.e. and in photons.

All its parameters are configurable per telescope, in the config files (lists, e.g.
``"gain_drop_slope": [["id", 21, -0.00487], ["id", 22, -0.00118]]``) or on the command line
(tuples, e.g. ``--SiPMVoltageDropModel.gain_drop_slope='[("id", 21, -0.00487), ("id", 22, -0.00118)]'``;
a single value applies to all the telescopes). As for every ctapipe telescope parameter, a
configured value replaces the whole default list: give the values of all the telescopes.

The functions read the statistics of the ADC samples of the pedestal events of monitoring or DL1
files and plot the NSB rate, see also the ``sst1mpipe-nsb`` tool.
"""
import logging

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import tables
from astropy.table import vstack
from ctapipe.core import TelescopeComponent
from ctapipe.core.traits import Dict, FloatTelescopeParameter, TelescopeParameter
from scipy import interpolate

__all__ = [
    "BASELINE_SHIFT_SPLINES",
    "SiPMVoltageDropModel",
    "mean_pedestal_std",
    "pedestal_time",
    "plot_nsb_rate",
    "read_pedestal_tables",
]

logger = logging.getLogger(__name__)

# B-splines (scipy representation) of the baseline shift (ADC) vs. the variance of the ADC
# samples of the pedestal events (ADC²), fitted for each telescope
BASELINE_SHIFT_SPLINES = {
    21: dict(
        knots=[0.0, 0.0, 0.0, 46.96969697, 60.90909091, 67.57575758, 74.24242424, 74.24242424, 74.24242424],
        coefficients=[-0.4367913, 9.08788226, 39.66927065, 54.18389959, 69.24518367, 104.29418015, 0.0, 0.0, 0.0],
        degree=2,
    ),
    22: dict(
        knots=[
            0.0, 0.0, 93.74735497, 126.46232902, 181.54405469, 256.80398702, 266.62290877,
            292.98812266, 294.71402901, 295.67072519, 298.04051867, 339.91032956, 339.91032956,
        ],
        coefficients=[
            -1.74764142, 48.41280122, 68.36843193, 104.62272377, 168.14318275, 179.34289593,
            212.07400565, 214.10474285, 215.17017626, 217.7828114, 268.19128994, 0.0, 0.0,
        ],
        degree=1,
    ),
}


class SiPMVoltageDropModel(TelescopeComponent):
    """
    Voltage drop of the SiPMs of the SST-1M cameras under NSB light, from the variance of the ADC
    samples of the pedestal events (ADC², the square of their std).

    - baseline shift (ADC): ``baseline_shift_spline`` of the variance;
    - gain drop: ``1 + gain_drop_slope * baseline shift``;
    - charge drop (gain and PDE): ``1 + (gain_drop_slope + pde_drop_slope) * baseline shift``;
    - NSB rate (p.e.): baseline shift / (``pulse_integral`` * gain drop) / duration of a sample;
    - NSB photon rate: baseline shift / (``pulse_integral`` * charge drop) / duration of a sample
      / ``photon_detection_efficiency``.

    The duration of a sample is given by the sampling rate of the camera (subarray). The telescopes
    without parameters (e.g. the simulated telescopes) raise a ValueError.
    """

    baseline_shift_spline = TelescopeParameter(
        trait=Dict(allow_none=True),
        default_value=[("type", "*", None)] + [("id", tel_id, spline) for tel_id, spline in BASELINE_SHIFT_SPLINES.items()],
        allow_none=True,
        help=(
            "B-spline of the baseline shift (ADC) vs. the variance of the ADC samples of the pedestal"
            " events (ADC²): dict of the knots, coefficients and degree (scipy splev). Out of the range"
            " of the knots, the spline is extrapolated."
        ),
    ).tag(config=True)

    gain_drop_slope = FloatTelescopeParameter(
        default_value=[("type", "*", None), ("id", 21, -0.00487), ("id", 22, -0.00118)],
        allow_none=True,
        help="Relative drop of the gain per ADC of baseline shift (arXiv:1910.00348, fig. 11.1)",
    ).tag(config=True)

    pde_drop_slope = FloatTelescopeParameter(
        default_value=[("type", "*", None), ("id", 21, -0.00317), ("id", 22, -0.000663)],
        allow_none=True,
        help="Relative drop of the PDE per ADC of baseline shift (arXiv:1910.00348, fig. 11.1)",
    ).tag(config=True)

    pulse_integral = FloatTelescopeParameter(
        default_value=[("type", "*", None), ("id", 21, 21.0), ("id", 22, 23.5)],
        allow_none=True,
        help="Integral of the pulse of one p.e. without voltage drop (ADC x samples)",
    ).tag(config=True)

    photon_detection_efficiency = FloatTelescopeParameter(
        default_value=0.35,
        help="PDE of the SiPMs without voltage drop, for the NSB photon rate (approximate)",
    ).tag(config=True)

    def _parameter(self, name, tel_id):
        value = getattr(self, name).tel[tel_id]
        if value is None:
            raise ValueError(f"No {type(self).__name__}.{name} for the telescope {tel_id}")
        return value

    def baseline_shift(self, tel_id, variance):
        """Baseline shift (ADC) for the variance of the ADC samples of the pedestal events (ADC²)"""
        spline = self._parameter("baseline_shift_spline", tel_id)
        tck = (np.asarray(spline["knots"], dtype=float), np.asarray(spline["coefficients"], dtype=float), int(spline["degree"]))
        return interpolate.splev(variance, tck)

    def gain_drop(self, tel_id, variance):
        """Factor by which the gain is reduced"""
        return 1 + self._parameter("gain_drop_slope", tel_id) * self.baseline_shift(tel_id, variance)

    def charge_drop(self, tel_id, variance):
        """Factor by which the charges are reduced (gain and PDE): the voltage drop correction divides the charges by it"""
        slope = self._parameter("gain_drop_slope", tel_id) + self._parameter("pde_drop_slope", tel_id)
        return 1 + slope * self.baseline_shift(tel_id, variance)

    def sample_duration(self, tel_id):
        """Duration of an ADC sample of the camera of the telescope"""
        return (1 / self.subarray.tel[tel_id].camera.readout.sampling_rate).to(u.ns)

    def nsb_rate(self, tel_id, variance):
        """NSB rate per pixel (p.e.)"""
        pulse_integral = self._parameter("pulse_integral", tel_id)
        shift = self.baseline_shift(tel_id, variance)
        rate = shift / (pulse_integral * self.gain_drop(tel_id, variance)) / self.sample_duration(tel_id)
        return rate.to(u.MHz)

    def nsb_photon_rate(self, tel_id, variance):
        """NSB photon rate per pixel"""
        pulse_integral = self._parameter("pulse_integral", tel_id)
        shift = self.baseline_shift(tel_id, variance)
        rate = shift / (pulse_integral * self.charge_drop(tel_id, variance)) / self.sample_duration(tel_id)
        return (rate / self.photon_detection_efficiency.tel[tel_id]).to(u.MHz)


def default_model(config=None):
    """SiPMVoltageDropModel of the SST-1M telescopes (subarray of sst1mpipe), configured by ``config``"""
    from sst1mpipe.utils.utils import get_subarray

    return SiPMVoltageDropModel(subarray=get_subarray(), config=config)


def read_pedestal_tables(files, last_chunk_only=False):
    """
    Statistics of the ADC samples of the pedestal events of monitoring or DL1 files (see
    `~sst1mpipe.io.load_r0_pedestals`, with the column tel_id), stacked in the order of the file
    names. ``last_chunk_only``: only the last chunk of each file and telescope. The files without
    these statistics are skipped (warning); None if no file has them.
    """
    from sst1mpipe.io import load_r0_pedestals

    pedestal_tables = []
    for path in sorted(map(str, files)):
        try:
            table = load_r0_pedestals(path)
        except tables.NoSuchNodeError:
            logger.warning("No statistics of the pedestal events in %s", path)
            continue
        if last_chunk_only:
            table = vstack([group[-1:] for group in table.group_by("tel_id").groups])
        pedestal_tables.append(table)
    return vstack(pedestal_tables) if pedestal_tables else None


def mean_pedestal_std(pedestal_table):
    """Mean over the camera pixels of the std of the ADC samples of the pedestal events (one per row)"""
    std = np.asarray(pedestal_table["std"])
    return np.nanmean(std.reshape(len(std), -1), axis=1)


def pedestal_time(pedestal_table):
    """Time of the statistics of the pedestal events: middle of the chunk of events (one per row)"""
    return pedestal_table["time_start"] + (pedestal_table["time_end"] - pedestal_table["time_start"]) / 2


def plot_nsb_rate(pedestal_table, tel_id, model=None, photon_rate=False, ax=None, **plot_kwargs):
    """
    NSB rate (p.e., or photons with ``photon_rate``) vs. time, from the variance of the ADC samples
    (mean over the pixels of their std, squared) of each row of ``pedestal_table`` (one telescope).
    ``plot_kwargs`` are given to matplotlib ``plot``. Returns the axes.
    """
    model = default_model() if model is None else model
    variance = mean_pedestal_std(pedestal_table) ** 2
    rate = model.nsb_photon_rate(tel_id, variance) if photon_rate else model.nsb_rate(tel_id, variance)
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 5))
    ax.plot(pedestal_time(pedestal_table).utc.to_datetime(), rate.to_value(u.MHz), ".", **plot_kwargs)
    ax.set_xlabel("Time")
    ax.set_ylabel(f"NSB {'photon rate' if photon_rate else 'rate'} [MHz]")
    return ax
