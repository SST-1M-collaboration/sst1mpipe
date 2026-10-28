.. _sst1mpipe_monitoring:

Monitoring of a run: ``sst1mpipe-monitoring``
=============================================

``sst1mpipe-monitoring`` is a ctapipe tool which computes the monitoring data of a run of a SST-1M
telescope from its raw data (R0, ``SST1M<tel>_<date>_<run>.fits.fz``): the pointing of the telescope and
the statistics of its pedestal events. The monitoring file is then given to ``sst1mpipe-process``
(``--monitoring-source SST1MMonitoringSource``), which uses it for the calibration and the cleaning of the
events of the run, see :ref:`data_processing`.

.. code-block:: bash

    sst1mpipe-monitoring \
        --input SST1M1_20260120_1179.fits.fz \
        --output SST1M1_20260120_1179.monitoring.h5 \
        --config sst1mpipe/resources/config/sst1mpipe_monitoring_config.json

Main options (``sst1mpipe-monitoring --help``, and ``--help-all`` for all the settings):

=============================  ===================================================================
``-i``, ``--input``            raw data file of the run (``EventSource.input_url``)
``-o``, ``--output``           monitoring file (HDF5), required
``--config``                   configuration file, by default
                               ``sst1mpipe/resources/config/sst1mpipe_monitoring_config.json``
``-m``, ``--max-events``       number of events read (all the events of the run by default)
``-t``, ``--allowed-tels``     telescopes read (all by default)
``--overwrite``                overwrite the output file (an existing file is an error otherwise)
``--progress``                 progress bar during the event loop
=============================  ===================================================================

Only the SST-1M raw data are accepted: the simulations (``.simtel.gz``) have no monitoring data and are
processed without monitoring source.


What the tool does
------------------

The raw data file is read twice, so that the images of the pedestal events are calibrated exactly as the
events by ``sst1mpipe-process``, with the statistics of the ADC samples:

1. **First pass:** the time of the first and of the last event of each telescope is kept, for the
   pointing table, and the mean and the std of the ADC samples of each pixel (R0 waveforms) of the
   **pedestal events** (internal triggers of the camera). The statistics of the ADC samples are then
   computed (rolling chunks, see below), and from them the **camera calibration coefficients**
   (``CameraCoefficientsCalculator``).
2. **Second pass:** the image of each pedestal event (charge of each pixel in p.e.) is calibrated as the
   events by ``sst1mpipe-process``:

   * the camera calibration coefficients of the first pass at the time of the pedestal event (last chunk
     ended before it, as given by ``SST1MMonitoringSource`` to ``sst1mpipe-process``);
   * ``R0R1Calibrator``: DigiCam pedestal, dc to p.e., window, bad pixels;
   * ``CameraCalibrator``: integration of the waveforms, multiplied by the ``factor`` of the coefficients
     (**voltage drop correction**), the charge of the **dead pixels** and of the bad pixels interpolated
     from their neighbours;
   * ``ImageSaturationCorrector``.

   The statistics of the images are then computed. The pedestal events without coefficients (no chunk of
   statistics of the ADC samples at their time) are not used (warning in the log), as they are skipped by ``sst1mpipe-process``.

The calibration settings (``R0R1Calibrator``, ``CameraCalibrator``, ``ImageSaturationCorrector``) of
``sst1mpipe_monitoring_config.json`` are those of ``sst1mpipe_data_config.json`` and
``sst1mpipe_rta_config.json``: if they are changed for the processing, they must be changed for the
monitoring too (checked by the tests). Reading the file twice takes about 25 % more time than once (43 s
instead of 35 s for 11.5 s of data of a telescope).

At the end of the run, the tool writes:

* the **pointing** of each telescope: altitude and azimuth of the target of the run every
  ``MonitoringTool.pointing_update_interval`` (1 s) from the first to the last event of the telescope.
  The target is read from the TARGET field of the file, or given with ``--SST1MEventSource.pointing_ra``
  and ``--SST1MEventSource.pointing_dec``. The runs without target (e.g. dark runs, transitions between
  two wobbles) have no pointing table, with a warning.
* the **statistics of the ADC samples of the pedestal events**, in rolling chunks of pedestal events
  (``WaveformStatisticsCalculator``, the last 100 pedestal events, at each pedestal event, by default): ``mean`` and ``median`` of the baseline (mean
  of the ADC samples), ``std`` the noise (mean of the std of the ADC samples, used for the voltage drop
  correction and the dead pixels), ``outlier_mask`` (dead pixels: noise below half of the median of the
  camera) and ``is_valid``.
* the **camera calibration coefficients**, one row per chunk of the statistics of the ADC samples, in the
  ctapipe table ``dl1/monitoring/telescope/calibration/camera/coefficients/tel_XXX`` (read by the ctapipe
  ``HDF5MonitoringSource``): ``factor`` = 1 / charge drop of the ``SiPMVoltageDropModel`` (voltage drop
  correction, ``CameraCoefficientsCalculator.voltage_drop_correction``: ``global`` with the median variance
  of the camera, ``pixelwise`` or ``none``), ``outlier_mask`` (dead pixels, ``flag_dead_pixels``),
  ``time`` (end of the chunk: valid from its last pedestal event), ``time_start``, ``time_end``,
  ``is_valid``, and ``pedestal_offset`` and ``time_shift`` (0).
* the **statistics of the images of the pedestal events**, in rolling chunks of pedestal events
  (``ImageStatisticsCalculator``, the last 1000 pedestal events, at each pedestal event, by default): ``mean``, ``median`` and ``std`` of the charges
  (p.e.). The ``std`` is the NSB noise used by the ``NSBImageCleaner`` to raise the picture threshold of
  the pixels, and the ``mean`` gives ``meanQ`` of the DL1 parameters.
* the description of the telescope, the scheduling and observation blocks of the run and the reference
  metadata, as in the ctapipe DL1 files.

The event type of the pedestal events, ``<event type>`` in the names of the tables, is given by the target
of the run:

=====================  =========================  ===================================
target of the run      event type                 tables
=====================  =========================  ===================================
``dark``               ``DARK_PEDESTAL``          ``dark_pedestal_waveform/image``
``bias``               ``ELECTRONIC_PEDESTAL``    ``electronic_pedestal_waveform/image``
other targets          ``SKY_PEDESTAL``           ``sky_pedestal_waveform/image``
=====================  =========================  ===================================

A run without pedestal events has no statistics tables (warning in the log), and a run without any event
gives no output file.


Rolling chunks of the last pedestal events
------------------------------------------

The statistics are computed in **rolling chunks**: a chunk of ``chunk_size`` pedestal events starts every
``chunk_shift`` pedestal events, so that consecutive chunks overlap. For each event, ``sst1mpipe-process``
uses the **last chunk ended before the event**: the statistics of the last ``chunk_size`` pedestal events
before it, updated every ``chunk_shift`` pedestal events. By default ``chunk_shift`` is 1: a new chunk at
each pedestal event, so that the statistics of an event are always those of the ``chunk_size`` pedestal
events just before it (100 for the ADC samples, 1000 for the images). For example, with ``chunk_size`` 4
and ``chunk_shift`` 1:

.. code-block:: text

    pedestal events   p0  p1  p2  p3  p4  p5  p6  p7
    chunk 0           [p0 ........... p3]
    chunk 1               [p1 ........... p4]
    chunk 2                   [p2 ........... p5]
    chunk 3                       [p3 ........... p6]
    chunk 4                           [p4 ........... p7]

    event between p3 and p4  ->  chunk 0  (p0-p3)
    event between p4 and p5  ->  chunk 1  (p1-p4)
    event between p6 and p7  ->  chunk 3  (p3-p6)
    event after p7           ->  chunk 4  (p4-p7)

With ``chunk_shift`` 2, there would be the chunks 0, 2 and 4 only: the events between p5 and p7 would
use chunk 2 (p2-p5).

* The events before the end of the first chunk (the first ``chunk_size`` pedestal events of the run) use
  the first chunk: there are not yet ``chunk_size`` pedestal events before them.
* The last chunk always ends at the last pedestal event of the run (``last_chunk_policy`` ``overlap``),
  and is not repeated if the last rolling chunk already ends there.
* The events before the first pedestal event or after the last one, by more than
  ``SST1MMonitoringSource.timestamp_tolerance`` (1 s), have no statistics and are skipped by
  ``sst1mpipe-process``.
* A run with less pedestal events than ``chunk_size`` has a single chunk of all its pedestal events
  (``allow_undersized_tables``).

A run of n pedestal events has n - ``chunk_size`` + 1 chunks: with about 100 pedestal events per second,
a run of 40 s has about 4000 pedestal events, about 3900 chunks of ADC samples (the last second of
pedestal events) and 3000 chunks of images (the last 10 s). A larger ``chunk_shift`` gives less rows in
the monitoring and DL1 files, but updates the statistics less often.

For each chunk, as in the ``ctapipe-calculate-pixel-statistics`` tool (ctapipe
``PixelStatisticsCalculator``, ``sst1mpipe.monitoring``):

* the statistics of each pixel are computed over the pedestal events of the chunk;
* the outlier pixels are flagged by the outlier detectors (``outlier_mask``), and are NaN in the
  statistics given to the events;
* the chunk is not valid (``is_valid``) if more than ``faulty_pixels_fraction`` (10 %) of its pixels are
  outliers.

Unlike the ``ctapipe-calculate-pixel-statistics`` tool, whose first pass ignores ``chunk_shift``, all the
chunks are rolling, and the chunk of an event is the last one ended before it, not the last one started
before it (which would contain the pedestal events after the event).

The main settings, in ``sst1mpipe_monitoring_config.json``:

* ``WaveformStatisticsCalculator.PedestalWaveformAggregator.SizeChunking`` and
  ``ImageStatisticsCalculator.PlainAggregator.SizeChunking``: ``chunk_size`` the number of pedestal
  events in a chunk, ``chunk_shift`` the number of pedestal events between the starts of two chunks (at
  most ``chunk_size``, ``null`` for chunks which do not overlap). ``chunking_type`` ``TimeChunking``
  gives chunks of a duration instead (``chunk_duration`` and ``chunk_shift`` in seconds).
  ``allow_undersized_tables`` gives a single chunk for the runs with less pedestal events than a chunk,
  and ``last_chunk_policy`` ``overlap`` makes the last chunk end at the last pedestal event.
* ``stats_aggregator_type`` of the ``ImageStatisticsCalculator``: ``PlainAggregator`` (mean, median and
  std) or ``SigmaClippingAggregator`` (without the outlier events, ``max_sigma`` and ``iterations``)
* ``outlier_detector_list``: ``RangeOutlierDetector`` (``validity_range``), ``MedianOutlierDetector``
  (``median_range_factors``) or ``StdOutlierDetector`` (``std_range_factors``) applied to the ``mean``,
  ``median`` or ``std`` of the chunks
* ``CameraCoefficientsCalculator``: ``voltage_drop_correction`` and ``flag_dead_pixels``, and the
  ``SiPMVoltageDropModel``, for the camera calibration coefficients
* ``R0R1Calibrator`` and ``CameraCalibrator``: calibration of the images, as for the processing
* ``SST1MEventSource.swapped_modules``: modules of the camera whose cables were swapped during a period,
  which are put back in place when the events are read

Every setting can also be given on the command line, e.g.
``--WaveformStatisticsCalculator.PedestalWaveformAggregator.SizeChunking.chunk_shift=10`` (a new chunk every
10 pedestal events).


Content of the monitoring file
------------------------------

.. code-block:: text

    SST1M1_20260120_1179.monitoring.h5
    ├── configuration/                    telescope description, scheduling and observation blocks
    ├── dl0/monitoring/telescope/pointing/tel_021
    │       time, azimuth, altitude
    ├── r0/monitoring/telescope/calibration/camera/pixel_statistics/sky_pedestal_waveform/tel_021
    │       time_start, time_end, event_id_start, event_id_end, n_events,
    │       mean, median, std, outlier_mask, is_valid
    └── dl1/monitoring/telescope/calibration/camera/pixel_statistics/sky_pedestal_image/tel_021
            time_start, time_end, event_id_start, event_id_end, n_events,
            mean, median, std, outlier_mask, is_valid

The monitoring file is a ctapipe monitoring file: it is read by the ctapipe ``HDF5MonitoringSource``
(pointing and statistics of the images), by ``SST1MMonitoringSource`` (all the tables) and its tables with
``ctapipe.io.read_table``:

.. code-block:: python

    from ctapipe.io import read_table

    path = "/dl1/monitoring/telescope/calibration/camera/pixel_statistics/sky_pedestal_image/tel_021"
    statistics = read_table("SST1M1_20260120_1179.monitoring.h5", path)
    print(statistics["time_start", "time_end", "n_events"])

The monitoring tables are also copied into the DL1 file by ``sst1mpipe-process``
(``--ProcessorTool.attach_monitoring``), so that the DL1 file is self-contained.
