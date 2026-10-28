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

The command line tools ``sst1mpipe-monitoring`` and ``sst1mpipe-process`` are installed with
``sst1mpipe`` (see :ref:`introduction`).


.. _processing_config:

Configuration
-------------

All the tools take a configuration file (``--config``, JSON or YAML), and every setting can also be given
on the command line, e.g. ``--WaveformStatisticsCalculator.PedestalWaveformAggregator.SizeChunking.chunk_size=50``.
``--help`` lists the main options and ``--help-all`` all the settings. The configuration used is stored
in the provenance log of the tool (``<tool>.provenance.log``, ``--provenance-log``).

The configuration files of ``sst1mpipe`` are in ``sst1mpipe/resources/config``:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - File
     - Used by
   * - ``sst1mpipe_monitoring_config.json``
     - ``sst1mpipe-monitoring``
   * - ``sst1mpipe_data_config.json``
     - ``sst1mpipe-process`` for the real data, and the analysis steps after DL1
   * - ``sst1mpipe_mc_config_low_nsb.json``, ``sst1mpipe_mc_config_high_nsb.json``
     - ``sst1mpipe-process`` for the simulations at low and high NSB, and the analysis steps after DL1
   * - ``sst1mpipe_rta_config.json``
     - ``sst1mpipe-process`` for the ZMQ streams of the real-time analysis

Each section of a file configures a component (a ctapipe class) by its name. The settings of a
telescope component can be given per telescope, as a list of rules: ``["type", "*", value]`` for all the
telescopes, ``["id", 21, value]`` for telescope 21. When several rules match a telescope, the last one
wins. A plain value applies to all the telescopes.

**Calibration** (``sst1mpipe_data_config.json``): ``null`` files are the default files of the telescope
(``sst1mpipe/resources/calibration``: averaged calibration parameters of 2023-2024, and the window
transmittance measured in the lab).

.. literalinclude:: ../sst1mpipe/resources/config/sst1mpipe_data_config.json
   :language: json
   :start-at: "R0R1Calibrator": {
   :end-before: "telescope_equivalent_focal_length"
   :dedent: 4

**Cleaning** (``sst1mpipe_data_config.json``): the ``NSBImageCleaner`` with the thresholds of each telescope.
The picture threshold of a pixel is raised to ``pedestal_factor`` times the std of its pedestal images
when this is above ``picture_threshold_pe``.

.. literalinclude:: ../sst1mpipe/resources/config/sst1mpipe_data_config.json
   :language: json
   :start-at: "ImageProcessor": {
   :end-before: "DataWriter"
   :dedent: 4

**Rolling chunks of the monitoring** (``sst1mpipe_monitoring_config.json``): statistics of the last 100
pedestal events, updated at each pedestal event (see :ref:`sst1mpipe_monitoring`).

.. literalinclude:: ../sst1mpipe/resources/config/sst1mpipe_monitoring_config.json
   :language: json
   :start-at: "WaveformStatisticsCalculator": {
   :end-before: "ImageStatisticsCalculator"
   :dedent: 2

**Overriding some settings**: ``--config`` can be given several times; each file overrides the
settings it gives, and a list (e.g. the rules of a telescope parameter) is replaced as a whole. E.g. the
calibration parameters of dark runs of the two telescopes, and chunks of 50 pedestal events, in
``my_config.yaml``:

.. code-block:: yaml

    R0R1Calibrator:
      calibration_file:
        - [id, 21, sst1mpipe/resources/calibration/dark_runs/20230917_0528_0529_Tel1_fitted_parameters.h5]
        - [id, 22, sst1mpipe/resources/calibration/dark_runs/20230921_0149_0150_Tel2_fitted_parameters.h5]
    WaveformStatisticsCalculator:
      PedestalWaveformAggregator:
        SizeChunking:
          chunk_size: 50

.. code-block:: bash

    sst1mpipe-monitoring --input SST1M1_20260120_1179.fits.fz --output SST1M1_20260120_1179.monitoring.h5 \
        --config sst1mpipe/resources/config/sst1mpipe_monitoring_config.json --config my_config.yaml

The full configuration files:

.. toggle::

    .. literalinclude:: ../sst1mpipe/resources/config/sst1mpipe_monitoring_config.json
       :language: json
       :caption: sst1mpipe_monitoring_config.json

.. toggle::

    .. literalinclude:: ../sst1mpipe/resources/config/sst1mpipe_data_config.json
       :language: json
       :caption: sst1mpipe_data_config.json


.. _monitoring_step:

1. Monitoring of the run: ``sst1mpipe-monitoring``
--------------------------------------------------

.. code-block:: bash

    sst1mpipe-monitoring \
        --input SST1M1_20260120_1179.fits.fz \
        --output SST1M1_20260120_1179.monitoring.h5 \
        --config sst1mpipe/resources/config/sst1mpipe_monitoring_config.json

The raw data file is read twice (the images of the pedestal events are calibrated as the events by
``sst1mpipe-process``, with the camera calibration coefficients of the first pass), and the monitoring file
contains the pointing of the telescope, the statistics of the ADC samples of the pedestal events (baseline
and noise), the camera calibration coefficients computed from them (voltage drop correction and dead
pixels) and the statistics of their images (NSB noise, for the ``NSBImageCleaner``), in rolling chunks of
the last pedestal events. What the tool does, its settings and the content of the monitoring file are described
in :ref:`sst1mpipe_monitoring`.


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
  (``event.monitoring.tel[tel_id].camera.pixel_statistics.pedestal_waveform``)
* the camera calibration coefficients (``event.monitoring.tel[tel_id].camera.coefficients``): the
  ``factor`` of the voltage drop correction and the ``outlier_mask`` of the dead pixels (the outlier pixels
  of the monitoring), applied by the ctapipe ``CameraCalibrator``
* the statistics of the images of the pedestal events (``pixel_statistics.pedestal_image``), used by the
  ``NSBImageCleaner``

The statistics of an event are those of the last chunk of pedestal events ended before it: with the
rolling chunks of ``sst1mpipe-monitoring``, the statistics of the last 100 (ADC samples) and 1000
(images) pedestal events before the event, updated at each pedestal event. The events before
the end of the first chunk use the first chunk. The events out of the time range of the chunks, by more
than ``--SST1MMonitoringSource.timestamp_tolerance`` (1 s), or out of the pointing table, have no
monitoring data and cannot be calibrated: they are skipped, and
their number is given in the log. The pedestal events are calibrated and cleaned but not written:
they are used to count the fraction of pedestal events which survive the cleaning.

Then the events are calibrated (``R0R1Calibrator``, ``CameraCalibrator``, ``ImageSaturationCorrector``),
cleaned and parametrized (``ImageProcessor``), and their shower geometry is reconstructed
(``ShowerProcessor``) if the DL2 is written (``DataWriter.write_dl2``, true in the data config).

The output file contains the DL1 (and DL2) of the events, ``/dl1/info`` (summary of the processing:
calibration files, target, number of triggered and pedestal events, number of pedestal events which
survive the cleaning, average fraction of pixels whose picture threshold is raised by the
``NSBImageCleaner``) and, with
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

The ``R0R1Calibrator`` does the corrections which do not need the pedestal events: DigiCam pedestal,
dc to p.e., window transmittance, and the pixels with bad calibration parameters and the ``bad_pixels``
added to the ``outlier_mask``.

**Without monitoring source**, the raw data are calibrated without the coefficients of the pedestal
events: no voltage drop correction, no dead pixels (only the ``bad_pixels`` and the pixels with bad
calibration parameters, warning in the log), and the picture thresholds of the ``NSBImageCleaner`` are not raised. The
output file has no pointing (the event source gives the RA/Dec of the run only), so it cannot be
reconstructed in DL2 by ``sst1mpipe_dl1_dl2``. This is used for the ZMQ streams of the real-time analysis
(``sst1mpipe_rta_config.json``); a monitoring file can be given to the stream as well.

**Simulations** (``.simtel.gz`` files) are processed without monitoring source, with the MC
configuration files (``sst1mpipe_mc_config_low_nsb.json``, ``sst1mpipe_mc_config_high_nsb.json``). How
each telescope is calibrated is set in the configuration: for the simulated telescopes 1 and 2, the R1
waveforms of the simulation are kept (``R0R1Calibrator.calibrate_r0`` false), the charges are corrected for
the PDE drop by the ``factor`` of the coefficients (``R0R1Calibrator.pde_drop_factor``), and the saturated
pixels are not corrected (``ImageSaturationCorrector`` ``integration_level`` and ``peak_time_level``
``null``):

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


.. _nsb_level:

NSB level and voltage drop: ``sst1mpipe-nsb``
---------------------------------------------

The NSB light shifts the baseline of the SiPMs and increases the variance of their ADC samples. The NSB
current causes a voltage drop on the bias resistor of the SiPMs, which reduces their gain and their PDE
linearly with the baseline shift (fig. 11.1 of `SiPM behaviour under continuous light
<https://arxiv.org/abs/1910.00348>`_). The ``SiPMVoltageDropModel`` (``sst1mpipe.utils.nsb_tools``) gives,
from the variance of the ADC samples of the pedestal events (ADC², square of their std):

* the **baseline shift** (ADC): spline of the variance, fitted for each telescope (``baseline_shift_spline``);
* the **gain drop** ``1 + gain_drop_slope * shift`` and the **charge drop**
  ``1 + (gain_drop_slope + pde_drop_slope) * shift``: the voltage drop correction (``factor`` of the camera
  calibration coefficients of ``sst1mpipe-monitoring``) divides the charges by the charge drop;
* the **NSB rate** per pixel: shift / (``pulse_integral`` * gain drop) / duration of a sample (4 ns), and
  the **NSB photon rate**: with the charge drop and divided by ``photon_detection_efficiency``.

Its parameters are configurable per telescope in the config files (section ``SiPMVoltageDropModel``,
e.g. ``"gain_drop_slope": [["id", 21, -0.00487], ["id", 22, -0.00118]]``) and on the command line of the
tools which use it (``sst1mpipe-process``, ``sst1mpipe-monitoring``, ``sst1mpipe-nsb``), with tuples:
``--SiPMVoltageDropModel.gain_drop_slope='[("id", 21, -0.00487), ("id", 22, -0.00118)]'`` (a single value
applies to all the telescopes). As for every ctapipe telescope parameter, a value replaces the whole
default list, so give the values of all the telescopes. The spline is extrapolated out of the range of its
knots (variance above 74 ADC² for tel 21, 340 ADC² for tel 22).

``sst1mpipe-nsb`` writes, for each chunk of pedestal events of monitoring or DL1 files, the variance, the
baseline shift, the gain and charge drops and the NSB rates (ECSV, FITS or HDF5 table), and plots the NSB
rate vs. time:

.. code-block:: bash

    sst1mpipe-nsb -i SST1M1_20260120_1179.monitoring.h5 -i SST1M2_20260120_1102.monitoring.h5 \
        -o nsb.ecsv --plot nsb.png

(one ``-i`` per input file).


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
