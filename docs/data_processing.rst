.. _data_processing:

Processing the data with the ctapipe tools
==========================================

The raw data (R0) of the SST-1M telescopes are processed with two
`ctapipe <https://github.com/cta-observatory/ctapipe>`_ tools, run one after the other on each
run (raw data file ``SST1M<tel>_<date>_<run>.fits.fz`` of a telescope):

1. ``sst1mpipe-monitoring``: monitoring data of the run (pointing and statistics of the pedestal
   events), see :ref:`monitoring_step`
2. ``sst1mpipe-process``: calibration and DL1 (and DL2 geometry) of the events of the run, with the
   monitoring data of the run, see :ref:`process_step`

The DL1 files of the runs can then be merged, see :ref:`merging_step`.

.. code-block:: text

    SST1M1_20260120_1179.fits.fz ──► sst1mpipe-monitoring ──► SST1M1_20260120_1179.monitoring.h5
              │                                                          │
              └──────────────────► sst1mpipe-process ◄───────────────────┘
                                          │
                                          ▼
                            SST1M1_20260120_1179_W1.dl1.h5  ──►  ctapipe-merge

**Installation**: ``sst1mpipe`` needs ``ctapipe`` 0.32 and Python 3.12, see ``environment.yml``:

.. code-block:: bash

    mamba env create -f environment.yml
    mamba activate sst1m-dev
    pip install -e .

All the tools take a configuration file (``--config``, JSON or YAML), and every setting can also be given
on the command line, e.g. ``--WaveformStatisticsCalculator.PedestalWaveformAggregator.SizeChunking.chunk_size=50``.
``--help`` lists the main options and ``--help-all`` all the settings. The configuration used is stored
in the provenance log of the tool (``<tool>.provenance.log``, ``--provenance-log``).


.. _monitoring_step:

1. Monitoring of the run: ``sst1mpipe-monitoring``
--------------------------------------------------

.. code-block:: bash

    sst1mpipe-monitoring \
        --input SST1M1_20260120_1179.fits.fz \
        --output SST1M1_20260120_1179.monitoring.h5 \
        --config sst1mpipe/resources/config/sst1mpipe_monitoring_config.json

The raw data file is read once, and the monitoring file contains:

* **pointing** of the telescope, ``/dl0/monitoring/telescope/pointing/tel_XXX``: altitude and azimuth
  of the target of the run (TARGET field of the file, or ``--SST1MEventSource.pointing_ra`` and
  ``--SST1MEventSource.pointing_dec``) every ``--MonitoringTool.pointing_update_interval`` (1 s) from the first to the
  last event of the run. The runs without pointing (e.g. dark runs) have no pointing table.
* **statistics of the ADC samples of the pedestal events**,
  ``/r0/monitoring/telescope/calibration/camera/pixel_statistics/<event type>_waveform/tel_XXX``,
  in chunks of pedestal events (100 by default): ``mean`` and ``median`` of the baseline (mean of the
  ADC samples), ``std`` the noise (std of the ADC samples, used for the voltage drop correction and
  the dead pixels), ``outlier_mask`` (dead pixels: noise below half the camera median) and ``is_valid``.
* **statistics of the images of the pedestal events**,
  ``/dl1/monitoring/telescope/calibration/camera/pixel_statistics/<event type>_image/tel_XXX``,
  in chunks of pedestal events (1000 by default): ``mean``, ``median`` and ``std`` of the charges (p.e.),
  used by the ``NSBImageCleaner``. The images are calibrated by the ``R0R1Calibrator`` (without voltage
  drop correction) and the ``CameraCalibrator``.
* the description of the telescope, the scheduling and observation blocks of the run and the reference
  metadata.

The pedestal events are the internal triggers of the camera. Their event type, ``<event type>`` in the
names of the tables, is given by the target of the run:

=====================  =========================  ===================================
target of the run      event type                 tables
=====================  =========================  ===================================
``dark``               ``DARK_PEDESTAL``          ``dark_pedestal_waveform/image``
``bias``               ``ELECTRONIC_PEDESTAL``    ``electronic_pedestal_waveform/image``
other targets          ``SKY_PEDESTAL``           ``sky_pedestal_waveform/image``
=====================  =========================  ===================================

The statistics are computed as in the ``ctapipe-calculate-pixel-statistics`` tool, by the
``WaveformStatisticsCalculator`` and the ``ImageStatisticsCalculator`` (ctapipe
``PixelStatisticsCalculator``): aggregation of the pedestal events in chunks, detection of the outlier
pixels, and validity of the chunks (``faulty_pixels_fraction``). The chunk of the time of an event is
the chunk from its first to its last pedestal event. Their main settings, in
``sst1mpipe_monitoring_config.json``:

* ``WaveformStatisticsCalculator.PedestalWaveformAggregator.SizeChunking.chunk_size`` and
  ``ImageStatisticsCalculator.PlainAggregator.SizeChunking.chunk_size``: number of pedestal events in a
  chunk. ``chunking_type`` ``TimeChunking`` gives chunks of a duration (``chunk_duration``) instead.
  ``allow_undersized_tables`` gives a single chunk for the runs with less pedestal events than a chunk.
* ``stats_aggregator_type`` of the ``ImageStatisticsCalculator``: ``PlainAggregator`` (mean, median and
  std) or ``SigmaClippingAggregator`` (without the outlier events)
* ``outlier_detector_list``: ``RangeOutlierDetector`` (``validity_range``), ``MedianOutlierDetector``
  (``median_range_factors``) or ``StdOutlierDetector`` (``std_range_factors``) applied to the ``mean``,
  ``median`` or ``std`` of the chunks
* ``R0R1Calibrator`` and ``CameraCalibrator``: calibration of the images, as for the processing

The monitoring file is a ctapipe monitoring file: it is read by the ctapipe ``HDF5MonitoringSource``
(pointing and statistics of the images) and its tables with ``ctapipe.io.read_table``:

.. code-block:: python

    from ctapipe.io import read_table

    path = "/dl1/monitoring/telescope/calibration/camera/pixel_statistics/sky_pedestal_image/tel_021"
    statistics = read_table("SST1M1_20260120_1179.monitoring.h5", path)
    print(statistics["time_start", "time_end", "n_events"])


.. _process_step:

2. Processing of the run: ``sst1mpipe-process``
-----------------------------------------------

.. code-block:: bash

    sst1mpipe-process \
        --input SST1M1_20260120_1179.fits.fz \
        --output SST1M1_20260120_1179.dl1.h5 \
        --config sst1mpipe/resources/config/sst1mpipe_data_config.json \
        --monitoring-source SST1MMonitoringSource \
        --SST1MMonitoringSource.input_files SST1M1_20260120_1179.monitoring.h5

As ``ctapipe-process --monitoring-source``, the monitoring source ``SST1MMonitoringSource`` fills the
monitoring containers of each event with the monitoring data at its time, before its calibration:

* the pointing of the telescope (``event.monitoring.tel[tel_id].pointing``), interpolated
* the statistics of the ADC samples of the pedestal events
  (``event.monitoring.tel[tel_id].camera.pixel_statistics.pedestal_waveform``), used by the
  ``R0R1Calibrator`` for the voltage drop correction and the dead pixels (the outlier pixels of the
  monitoring are dead pixels)
* the statistics of the images of the pedestal events (``pixel_statistics.pedestal_image``), used by the
  ``NSBImageCleaner``

The statistics of an event are those of the chunk of pedestal events at its time, or at most
``--SST1MMonitoringSource.timestamp_tolerance`` (1 s) before or after it. The events without monitoring
data (no chunk at their time, or out of the pointing table) cannot be calibrated: they are skipped, and
their number is given in the log. The pedestal events are not processed (they are used by the
monitoring only).

Then the events are calibrated (``R0R1Calibrator``, ``CameraCalibrator``, ``ImageSaturationCorrector``),
cleaned and parametrized (``ImageProcessor``), and their shower geometry is reconstructed
(``ShowerProcessor``) if the DL2 is written (``DataWriter.write_dl2``, true in the data config).

The output file contains the DL1 (and DL2) of the events, ``/dl1/info`` (summary of the processing:
calibration files, target, number of triggered and pedestal events) and, with
``--ProcessorTool.attach_monitoring=True`` (default), a copy of the monitoring data of the monitoring
file (``attach_monitoring_file``). The DL1 file is then self-contained: e.g. ``load_dl1_sst1m`` reads the pointing of the events
from it. The wobble of the run is added to the name of the output file (``_W1``,
``--ProcessorTool.wobble_in_output_name``).

Other settings of the processing:

* ``--ProcessorTool.allowed_sb_types``: the runs processed, by type of scheduling block
  (``OBSERVATION`` by default; ``UNKNOWN`` for the transitions between two wobbles, ``CALIBRATION`` for
  the dark and bias runs)
* ``--max-events``, ``--allowed-tels``: events and telescopes processed
* ``--DataWriter.write_dl1_images``, ``--DataWriter.write_r1_waveforms``: images and waveforms written
* ``--progress``: progress bar

**Without monitoring source**, the raw data are calibrated without the statistics of the pedestal
events: no voltage drop correction, no dead pixels (only the ``bad_pixels`` and the pixels with bad
calibration parameters), and the picture thresholds of the ``NSBImageCleaner`` are not raised. The
output file has no pointing (the event source gives the RA/Dec of the run only), so it cannot be
reconstructed in DL2 by ``sst1mpipe_dl1_dl2``. This is used for the ZMQ streams of the real-time analysis
(``sst1mpipe_rta_config.json``); a monitoring file can be given to the stream as well.

**Simulations** (``.simtel.gz`` files) are processed without monitoring source, with the MC
configuration files (``sst1mpipe_mc_config_low_nsb.json``, ``sst1mpipe_mc_config_high_nsb.json``):

.. code-block:: bash

    sst1mpipe-process --input gamma.simtel.gz --output gamma.dl1.h5 \
        --config sst1mpipe/resources/config/sst1mpipe_mc_config_low_nsb.json

**A night of data**: the monitoring and the processing of each run, e.g.

.. code-block:: bash

    CONFIG=sst1mpipe/resources/config
    for raw in /data/raw/2026/01/20/SST1M1/SST1M1_20260120_*.fits.fz; do
        run=$(basename "$raw" .fits.fz)
        sst1mpipe-monitoring -i "$raw" -o "$run.monitoring.h5" --config $CONFIG/sst1mpipe_monitoring_config.json
        sst1mpipe-process -i "$raw" -o "$run.dl1.h5" --config $CONFIG/sst1mpipe_data_config.json \
            --monitoring-source SST1MMonitoringSource --SST1MMonitoringSource.input_files "$run.monitoring.h5"
    done

The transition runs between two wobbles are skipped by ``sst1mpipe-process`` (no output events), and the
dark runs have no pointing.


.. _merging_step:

3. Merging
----------

**Monitoring data into a DL1 file**: if the DL1 file was produced without the monitoring data
(``--ProcessorTool.attach_monitoring=False``), they can be attached to it afterwards. The ctapipe
monitoring data (pointing, statistics of the images) are attached by ``ctapipe-merge`` with the
``monitoring-only`` merge strategy, and ``attach_monitoring_file`` also attaches the statistics of the
ADC samples (``/r0/monitoring``):

.. code-block:: bash

    ctapipe-merge --append --merge-strategy monitoring-only \
        SST1M1_20260120_1179.monitoring.h5 --output SST1M1_20260120_1179_W1.dl1.h5

.. code-block:: python

    from sst1mpipe.io import attach_monitoring_file

    attach_monitoring_file("SST1M1_20260120_1179.monitoring.h5", "SST1M1_20260120_1179_W1.dl1.h5")

**DL1 files of several runs of a telescope** (e.g. a night) are merged with ``ctapipe-merge`` (the files of the two telescopes have different subarrays: they are combined by the stereo step):

.. code-block:: bash

    ctapipe-merge --input-dir dl1/ --pattern "SST1M1_20260120_*.dl1.h5" --output SST1M1_20260120.dl1.h5
    # or
    ctapipe-merge SST1M1_20260120_1179_W1.dl1.h5 SST1M1_20260120_1180_W1.dl1.h5 --output merged.dl1.h5

The events, the pointing tables, the statistics of the images of the pedestal events and the
scheduling and observation blocks of the runs are merged. The statistics of the ADC samples
(``/r0/monitoring``) and ``/dl1/info`` are not: they stay in the monitoring and DL1 files of the runs.
``--no-dl1-images`` merges the files without the images. The simulated DL1 files are merged in the same
way (``sst1mpipe_merge_hdf5`` for a directory of MC files).


Reading the outputs
-------------------

.. code-block:: python

    from ctapipe.io import HDF5EventSource, TableLoader
    from sst1mpipe.io import load_dl1_sst1m, load_dl1_pedestals, load_r0_pedestals

    # DL1 parameters with the pointing and the trigger time of the events
    events = load_dl1_sst1m("SST1M1_20260120_1179_W1.dl1.h5", tel="tel_021")

    # statistics of the pedestal events of the run
    images = load_dl1_pedestals("SST1M1_20260120_1179_W1.dl1.h5", tel="tel_021")
    samples = load_r0_pedestals("SST1M1_20260120_1179.monitoring.h5", tel="tel_021")

    # ctapipe readers: the pointing of the events is interpolated from the pointing table
    with TableLoader("SST1M1_20260120_1179_W1.dl1.h5") as loader:
        table = loader.read_telescope_events([21], dl1_parameters=True, pointing=True)

``load_dl1_sst1m`` adds the parameters of the analysis which are not computed by ctapipe
(``add_dl1_extra_parameters``): ``log_camera_frame_hillas_intensity``, ``camera_frame_hillas_wl``,
``equivalent_focal_length``; for the observed data ``local_time``, ``date`` and ``meanQ`` (mean charge
of the pedestal images); for the simulations the true parameters of the showers, ``log_true_energy``,
the DISP parameters, ``miss`` and ``true_tel_impact_distance``. They are not written in the DL1 files.
