"""
Camera calibration coefficients of the SST-1M telescopes, from the statistics of the ADC samples of
the pedestal events: the corrections which need the pedestal events (voltage drop, dead pixels), as
ctapipe `~ctapipe.containers.CameraCalibrationContainer` coefficients applied by the ctapipe
`~ctapipe.calib.CameraCalibrator`. The corrections which do not need them are done by the
`~sst1mpipe.calib.R0R1Calibrator`.
"""
import numpy as np
from astropy.table import QTable
from ctapipe.core import TelescopeComponent
from ctapipe.core.traits import BoolTelescopeParameter, CaselessStrEnum, TelescopeParameter

from sst1mpipe.utils.nsb_tools import SiPMVoltageDropModel

__all__ = ["CameraCoefficientsCalculator", "VOLTAGE_DROP_CORRECTIONS"]

VOLTAGE_DROP_CORRECTIONS = ("none", "global", "pixelwise")


class CameraCoefficientsCalculator(TelescopeComponent):
    """
    Camera calibration coefficients of each chunk of the statistics of the ADC samples of the
    pedestal events (sst1mpipe-monitoring), written in the monitoring file and filled in
    ``event.monitoring.tel[tel_id].camera.coefficients`` by the
    `~sst1mpipe.io.monitoring_source.SST1MMonitoringSource`:

    - ``factor``: 1 / charge drop due to the voltage drop of the SiPMs with the NSB
      (``voltage_drop_correction``, `~sst1mpipe.utils.nsb_tools.SiPMVoltageDropModel`), from the std
      of the ADC samples: ``global`` with the median variance of the camera, ``pixelwise`` with the
      variance of each pixel, ``none`` 1. The ctapipe CameraCalibrator multiplies the charges by it.
    - ``outlier_mask``: the dead pixels (``flag_dead_pixels``), i.e. the pixels without statistics
      (outliers of the monitoring), whose charge is interpolated by the CameraCalibrator
    - ``pedestal_offset`` and ``time_shift``: 0 (the pedestal is subtracted by the R0R1Calibrator)
    """

    voltage_drop_correction = TelescopeParameter(
        trait=CaselessStrEnum(VOLTAGE_DROP_CORRECTIONS),
        default_value="global",
        help=(
            "Correction of the voltage drop due to the NSB: 'global' uses the median"
            " of the pedestal variance of the camera, 'pixelwise' the variance of each pixel"
        ),
    ).tag(config=True)

    flag_dead_pixels = BoolTelescopeParameter(
        default_value=True,
        help=(
            "Flag the dead pixels, whose charge is interpolated: the pixels without statistics of"
            " the ADC samples of the pedestal events (NaN), i.e. the outliers of the monitoring"
        ),
    ).tag(config=True)

    def __init__(self, subarray, config=None, parent=None, **kwargs):
        super().__init__(subarray=subarray, config=config, parent=parent, **kwargs)
        self.voltage_drop_model = SiPMVoltageDropModel(subarray=subarray, parent=self)

    def charge_drop(self, tel_id, pedestal_std):
        """
        Charge drop (factor by which the charges are reduced) of the pixels, from the std of the ADC
        samples of the pedestal events, shape (..., n_channels, n_pixels), NaN for the outliers
        """
        correction = self.voltage_drop_correction.tel[tel_id]
        variance = np.asarray(pedestal_std, dtype=np.float64) ** 2
        if correction == "none":
            return np.ones_like(variance)
        if correction == "global":
            # without the pixels without statistics (outliers of the monitoring)
            median = np.nanmedian(variance, axis=(-2, -1), keepdims=True)
            return np.broadcast_to(self.voltage_drop_model.charge_drop(tel_id, median), variance.shape).copy()
        return self.voltage_drop_model.charge_drop(tel_id, variance)

    def __call__(self, statistics, tel_id):
        """
        Coefficients of the chunks of ``statistics``, the statistics of the ADC samples of the
        pedestal events of the telescope (`~sst1mpipe.monitoring.compute_pixel_statistics`).

        Returns
        -------
        astropy.table.QTable
            one row per chunk, the columns of the ctapipe coefficients table (time, pedestal_offset,
            factor, time_shift, outlier_mask, is_valid) and the start and end of the chunks
            (time_start, time_end). The time is the end of the chunk: the coefficients are valid
            from the last pedestal event of the chunk, as the statistics (last chunk ended before
            an event).
        """
        std = np.array(statistics["std"], dtype=np.float64)
        std[np.asarray(statistics["outlier_mask"], dtype=bool)] = np.nan
        dead = np.isnan(std)
        if not self.flag_dead_pixels.tel[tel_id]:
            dead = np.zeros_like(dead)
        zeros = np.zeros(std.shape, dtype=np.float32)
        return QTable(dict(
            time=statistics["time_end"],
            time_start=statistics["time_start"],
            time_end=statistics["time_end"],
            pedestal_offset=zeros,
            factor=(1 / self.charge_drop(tel_id, std)).astype(np.float32),
            time_shift=zeros.copy(),
            outlier_mask=dead,
            is_valid=np.asarray(statistics["is_valid"], dtype=bool),
        ))
