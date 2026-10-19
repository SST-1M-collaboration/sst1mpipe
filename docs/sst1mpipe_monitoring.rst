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

The raw data file is read once, event by event:

1. the time of the first and of the last event of each telescope is kept, for the pointing table;
2. the **pedestal events** (internal triggers of the camera) are kept, the other events are only counted.
   For each pedestal event:

   * the mean and the std of the ADC samples of each pixel (R0 waveforms);
   * the image of the event (charge of each pixel in p.e.), calibrated as in ``sst1mpipe-process``
     by the ``R0R1Calibrator`` and the ``CameraCalibrator``, but without voltage drop correction and without
     dead pixels, since these need the statistics of the ADC samples which are being computed.

At the end of the run, the tool writes:

* the **pointing** of each telescope: altitude and azimuth of the target of the run every
  ``MonitoringTool.pointing_update_interval`` (1 s) from the first to the last event of the telescope.
  The target is read from the TARGET field of the file, or given with ``--SST1MEventSource.pointing_ra``
  and ``--SST1MEventSource.pointing_dec``. The runs without target (e.g. dark runs, transitions between
  two wobbles) have no pointing table, with a warning.
* the **statistics of the ADC samples of the pedestal events**, in chunks of pedestal events
  (``WaveformStatisticsCalculator``, 100 events by default): ``mean`` and ``median`` of the baseline (mean
  of the ADC samples), ``std`` the noise (mean of the std of the ADC samples, used for the voltage drop
  correction and the dead pixels), ``outlier_mask`` (dead pixels: noise below half of the median of the
  camera) and ``is_valid``.
* the **statistics of the images of the pedestal events**, in chunks of pedestal events
  (``ImageStatisticsCalculator``, 1000 events by default): ``mean``, ``median`` and ``std`` of the charges
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


Statistics in chunks of pedestal events
---------------------------------------

The statistics are computed as in the ``ctapipe-calculate-pixel-statistics`` tool, by the ctapipe
``PixelStatisticsCalculator`` (``sst1mpipe.monitoring``):

* the pedestal events are aggregated in chunks, and the statistics of each pixel are computed for each
  chunk. A chunk covers the time from its first to its last pedestal event: ``sst1mpipe-process`` uses
  for each event the chunk containing its time (within ``SST1MMonitoringSource.timestamp_tolerance``);
* the outlier pixels of each chunk are flagged by the outlier detectors (``outlier_mask``);
* a chunk is not valid (``is_valid``) if more than ``faulty_pixels_fraction`` (10 %) of its pixels are
  outliers. With overlapping chunks (``chunk_shift``), the time ranges of the invalid chunks are computed
  again with the overlapping chunks (second pass of ctapipe).

The main settings, in ``sst1mpipe_monitoring_config.json``:

* ``WaveformStatisticsCalculator.PedestalWaveformAggregator.SizeChunking.chunk_size`` and
  ``ImageStatisticsCalculator.PlainAggregator.SizeChunking.chunk_size``: number of pedestal events in a
  chunk. ``chunking_type`` ``TimeChunking`` gives chunks of a duration (``chunk_duration``) instead.
  ``allow_undersized_tables`` gives a single chunk for the runs with less pedestal events than a chunk,
  and ``last_chunk_policy`` ``overlap`` makes the last chunk overlap the previous one instead of being
  smaller.
* ``stats_aggregator_type`` of the ``ImageStatisticsCalculator``: ``PlainAggregator`` (mean, median and
  std) or ``SigmaClippingAggregator`` (without the outlier events, ``max_sigma`` and ``iterations``)
* ``outlier_detector_list``: ``RangeOutlierDetector`` (``validity_range``), ``MedianOutlierDetector``
  (``median_range_factors``) or ``StdOutlierDetector`` (``std_range_factors``) applied to the ``mean``,
  ``median`` or ``std`` of the chunks
* ``R0R1Calibrator`` and ``CameraCalibrator``: calibration of the images, as for the processing
* ``SST1MEventSource.swapped_modules``: modules of the camera whose cables were swapped during a period,
  which are put back in place when the events are read

Every setting can also be given on the command line, e.g.
``--WaveformStatisticsCalculator.PedestalWaveformAggregator.SizeChunking.chunk_size=50``.


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
