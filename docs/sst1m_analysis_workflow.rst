.. _sst1m_analysis_workflow:

SST-1M analysis workflow
========================

.. _automatic_processing:

Automatic data processing
-------------------------

Raw data from both telescopes are processed automaticaly up to DL3 level 
every morning after an observing night on **Calculus** with ``Daily analysis``
pipeline:

``/mnt/nfs22_auger/jurysek/sst1mpipe/sst1mpipe/scripts/daily_analysis.py``

Products of automatic processing (all data levels up to DL3) can be found in:

``/data/work/analysis/Daily_analysis/daily_data/``

Those are sorted by date, source observed, and the telescope used (``cs1``, ``cs2`` and ``stereo`` for coincident events).
All data level directories also contain configuration file and log file for each processed run, so 
user can easily check if the processing went OK and which version of ``sst1mpipe`` was used. Random Forest models used in authomatic analysis 
can be found in

``/data/work/analysis/MC/prod_january_2023/``

for both mono and stereo regimes of observation, together with expected performances and sensitivity curves. IRFs used for DL3
production can be found in

``/data/work/analysis/IRFs/``

DL3 products of the automatic processing can be further analyzed using standard tools 
such as `gammapy <https://github.com/gammapy>`_, see :ref:`higher_level`.

.. note::

    As the software development is progressing very fast, automatic processing is usually 
    done with a slightly obsolete version of ``sst1mpipe`` and with a configuration not fully 
    tunned for atmospheric conditions during particular nights. To obtain *paper level* 
    analysis, one should proceed with the steps described in :ref:`detailed_processing`.


.. _detailed_processing:

Detailed data processing
------------------------

All data files are stored on telescope servers cs1/cs2. They are also mirrored on Calculus
(``/net/cs1/data/raw/``, ``/net/cs2/data/raw/``) and `FZU farm <https://www.farm.particle.cz/>`_ (``/mnt/sst1m/sst1m_backups/``), 
where can be processed and analyzed.

To get some feeling of how the event reconstruction and data analysis works, one can take 
a look at materials from the last SST-1M hands-on analysis session in Prague, December 2023:

* `Introductory talk <https://indico.cern.ch/event/1337334/contributions/5692346/attachments/2775295/4836434/data_analysis_basics.pdf>`_
* Jupyter notebooks on Calculus: ``/data/work/analysis/sst1m_analysis_handson_2023/``

To understand general logic of the data processing, see scheme of the pipeline for mono and stereo
data analysis: :ref:`analysis_basics`.

Configuration files
~~~~~~~~~~~~~~~~~~~

There are two slightly different configuration files. One is to be used for MC processing
(needed only if one wants to get custom RF models / IRFs / performance curves), and one 
for real data processing. These contain configuration for all reconstruction/analysis steps 
from R0 to DL3 (or performance evaluation in case of MC).

Config files for MC, for the simulations at low and high NSB (they differ by the PDE drop correction factors):

.. toggle:: 

    .. include:: ../sst1mpipe/resources/config/sst1mpipe_mc_config_low_nsb.json
       :code: json

.. toggle:: 

    .. include:: ../sst1mpipe/resources/config/sst1mpipe_mc_config_high_nsb.json
       :code: json

Default config file for data:

.. toggle:: 

    .. include:: ../sst1mpipe/resources/config/sst1mpipe_data_config.json
       :code: json


R0 to DL1
~~~~~~~~~

The raw data (R0) are monitored with the ctapipe tool ``sst1mpipe-monitoring`` (statistics of the
pedestal events and pointing of each run, see :ref:`sst1mpipe_monitoring`), then the raw data (R0) and the
simulations (R1) are processed to DL1 with ``sst1mpipe-process``, see :ref:`data_processing`. The output is a ``HDF5`` file with a table of DL1
parameters, one output file per run (about 40 seconds of data) and per telescope (DL1 files from individual
telescopes searched for coincidences and merged for stereoscopic reconstruction in the next step
:ref:`dl1_dl1_stereo`).

It applies dc to p.e. calibration on raw R0 waveforms, integrates them, cleans the images (gets rid of the noisy pixels which likely do not
contain any Cherenkov photons), and parametrizes the shower images with Hillas ellipses.

**Relevant parts of the config file** applied in this analysis step:

* ``R0R1Calibrator`` - R0 -> R1 calibration (``sst1mpipe.calib.R0R1Calibrator``): the corrections which do not need the pedestal events. All settings can be given per telescope, e.g. ``[["type", "*", []], ["id", 22, [3, 1201]]]``:

  * ``calibration_file`` - calibration files based on analysis of dark runs (``null``: default file of the telescope). Should be taken relatively close to the date of observation
  * ``window_transmittance_file`` - camera window transmittance correction factors of the pixels, measured in the lab (``null``: default file of the telescope, can be kept)
  * ``flag_bad_calibration_pixels`` - the pixels with bad calibration parameters are flagged as outliers and interpolated
  * ``bad_pixels`` - ids of pixels always flagged as outliers and interpolated (data and MC)
  * ``pde_drop_factor`` - correction of the PDE drop due to the NSB in the simulations: the charges are divided by this factor (``factor`` of the camera calibration coefficients), which must match the PDE file of the simulation (``mc_pde_correction_factors.json``). ``null`` (real telescopes 21 and 22) applies no correction. ``sst1mpipe_mc_config_low_nsb.json`` and ``sst1mpipe_mc_config_high_nsb.json`` have the factors of the low and high NSB simulations

  The corrections which need the pedestal events are the camera calibration coefficients computed by ``sst1mpipe-monitoring`` (``CameraCoefficientsCalculator`` of the monitoring config, see :ref:`sst1mpipe_monitoring`): ``voltage_drop_correction`` (``none``, ``global`` or ``pixelwise`` correction of the voltage drop due to the NSB: the charges are divided by the charge drop of the ``SiPMVoltageDropModel``, see :ref:`nsb_level`) and ``flag_dead_pixels`` (the pixels without statistics of the ADC samples of the pedestal events are flagged as outliers and interpolated).

  The former ``telescope_calibration``, ``NsbCalibrator.*``, ``analysis.bad_pixels`` and ``window_transmittance`` settings are no longer read.

* ``NsbCalibrator.intensity_correction`` - global scale of the Hillas intensity per telescope. It is not applied at this step, but by ``load_dl1_sst1m(scale_intensities=True)``; not yet by ``sst1mpipe-apply-models`` (TODO)


* ``CameraCalibrator`` - Pulse integration settings

* ``ImageProcessor`` - Settings of image cleaning method and tailcuts. Tailcuts can be set different for different levels of NSB

* ``ShowerProcessor`` - Method of shower geometry reconstruction. Only applied if event source contains data from more telescopes, i.e. it's only relevant for MC in this analysis step.


Model training
~~~~~~~~~~~~~~

The random forests (energy regressor, particle classifier, disp reconstructor) are trained with the
ctapipe training tools (``sst1mpipe-train-energy-regressor``, ``sst1mpipe-train-particle-classifier``,
``sst1mpipe-train-disp-reconstructor``) on DL1 diffuse gammas and diffuse protons, see
:doc:`ml_reconstruction`.


.. _dl1_dl1_stereo:

DL1 to DL1 stereo
~~~~~~~~~~~~~~~~~

The telescopes have no array trigger: each telescope writes its own events, one raw data (zfits) file per
run and per telescope, processed into one DL1 file each by ``sst1mpipe-process``. The DL1 files of the runs
of a night are merged for each telescope, and the coincident events of the telescopes are combined into
stereo events by ``sst1mpipe-stereo``, a `ctapipe.core.Tool`:

.. code-block:: bash

    sst1mpipe-stereo \
        --input SST1M1_20260120_merged.dl1.h5 --input SST1M2_20260120_merged.dl1.h5 \
        --output SST1M_20260120_merged.stereo.dl1.h5

The matching method is chosen with ``--combiner`` (``--help`` lists the options). See
:ref:`stereo_combination`.

.. note::

    ``sst1mpipe-stereo`` **does not run on MC**: it stops with an error if an input file is simulated.


.. _dl1_dl2:

DL1 to DL2
~~~~~~~~~~

``sst1mpipe-apply-models`` applies the models of an NSB level (pointing of the telescopes as
feature) to the DL1 file (mono DL1 of ``sst1mpipe-process``, or stereo DL1 of
``sst1mpipe-stereo``, see :ref:`dl1_dl1_stereo`) and writes the ctapipe DL2, see :doc:`ml_reconstruction`.

.. note::

    The DL2 -> DL3 and IRF steps below read the DL2 of the former sst1mpipe random forests; they are
    to be migrated to the ctapipe DL2.


.. _dl2_irfs:

DL2 MC to IRFs
~~~~~~~~~~~~~~

.. note::

    In most cases, **analyser does not need to produce own Instrument Response Functions and this step can be safely skipped** using default IRFs referenced in 
    :ref:`automatic_processing`. IRF production, however, is necessary in performance studies, or if one uses custom RFs to produce DL2, or applies
    custom selection cuts in DL2 to DL3 step.

The IRFs are made with the tools of ctapipe from the DL2 of ``sst1mpipe-apply-models`` (point gammas,
diffuse gammas for the full enclosure IRFs, diffuse protons):

1. ``ctapipe-optimize-event-selection``: gammaness and direction (theta) cuts per reconstructed energy
   bin (``sst1mpipe_optimize_cuts.yaml``: the cuts of the best sensitivity; with
   ``sst1mpipe_irf_paper_cuts.yaml``: the fixed efficiencies of the SST-1M Crab paper, 60 % of the
   gammas for the gammaness cut, 80 % for the theta cut; with ``sst1mpipe_irf_strict_stereo.yaml``:
   events with both telescopes only), written to ``cuts.fits``
2. ``ctapipe-compute-irf`` (``sst1mpipe_compute_irf.yaml``): the IRFs with these cuts (effective area,
   energy dispersion, PSF or ``RAD_MAX``, background) in the GADF format of gammapy, and the
   benchmarks (energy bias and resolution, angular resolution, sensitivity in 50 h).
   Point-like IRFs: ``--spatial-selection-applied --DefaultFoVOffsetBins.fov_offset_n_bins=1``.

.. code-block:: console

    ctapipe-optimize-event-selection --config sst1mpipe_optimize_cuts.yaml \
        -g gamma_point.dl2.h5 -p proton_test.dl2.h5 -o cuts.fits
    ctapipe-compute-irf --config sst1mpipe_compute_irf.yaml --cuts cuts.fits \
        --gamma-file gamma_point.dl2.h5 --proton-file proton_test.dl2.h5 --do-background \
        --spatial-selection-applied --DefaultFoVOffsetBins.fov_offset_n_bins=1 \
        --output irf_point_like.fits.gz --benchmark-output benchmarks_point_like.fits.gz

The IRF files are compatible with gammapy:

.. code-block:: python

    from gammapy.irf import EffectiveAreaTable2D, EnergyDispersion2D

    aeff = EffectiveAreaTable2D.read("irf_point_like.fits.gz", hdu="EFFECTIVE AREA")
    edisp = EnergyDispersion2D.read("irf_point_like.fits.gz", hdu="ENERGY DISPERSION")

The figures of the IRFs and of the benchmarks are made by ``sst1mpipe-plot-irfs``, see :ref:`performance`.


.. _dl2_dl3:

DL2 to DL3
~~~~~~~~~~

``sst1mpipe_data_dl2_dl3`` is a tool to create DL3 data files from DL2 data files. It is supposed to be provided with a directory with input DL2 files 
(typicaly a directory with DL2 for one source observed in one night, but can be run on larger sample as well). It merges the DL2 ``HDF5`` per-run 
files into per-wobble DL3 ``fits`` files containing only photon lists. It also finds proper IRF based on zenith, azimuth and NSB level for each input DL2 file. It creates 
per-night index files needed for further analysis in gammapy. The script applies event selection cuts defined in the config file (``event_selection``), 
including cut on gammaness, so **one should make sure that these are the same as used on MC to produce IRFs**. Energy dependent gammaness cuts can be used as well, following 
the same rules described in :ref:`dl2_irfs`.

.. warning::

    ``sst1mpipe_data_dl2_dl3`` (and ``sst1mpipe_data_get_photon_list``) read the DL2 of the former
    random forests of sst1mpipe, not the DL2 of ctapipe written by ``sst1mpipe-apply-models``: the
    DL3 of the ctapipe DL2 is not implemented yet.


.. _higher_level:

High level analysis
~~~~~~~~~~~~~~~~~~~

Output DL3 files produced with ``sst1mpipe_data_dl2_dl3`` are fully compatible with gammapy and may be further analyzed using gammapy tools. See e.g. 

* `Tutorial on 1D spectral analysis <https://docs.gammapy.org/1.2/tutorials/analysis-1d/spectral_analysis.html>`_

* `Tutorial on 2D ring background map <https://docs.gammapy.org/1.2/tutorials/analysis-2d/ring_background.html>`_

A typical use case is to run joint gammapy analysis on data from several nights. In such case one has to run ``create_hdu_indexes`` script to create 
HDU index files indexing all DL3s to be used in the final analysis.


.. _performance:

Performance and sensitivity
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The performance is in the IRFs and benchmarks of ``ctapipe-compute-irf`` (see :ref:`dl2_irfs`), and
the performance of the models in the figures of the training tools, see :doc:`ml_reconstruction`.
``sst1mpipe-plot-irfs`` makes one figure per quantity (effective area, energy bias and resolution,
energy migration, angular resolution, cuts, background rate, sensitivity in 50 h, and with the DL2
the gammaness distributions and the ROC curve), with the same axis limits for all the figures.
Each ``--irf <label>=<directory>`` (``irf_point_like.fits.gz``, ``benchmarks_point_like.fits.gz``,
``cuts.fits``) is a series, e.g. to compare the mono analysis of each telescope and the stereo:

.. code-block:: console

    sst1mpipe-plot-irfs --irf tel1=irf/tel1 --irf tel2=irf/tel2 --irf stereo=irf/stereo \
        --dl2 stereo=dl2/stereo -o plots

``--reference <label>=<directory>`` draws IRF directories of the same format as lines: the
published performance is in ``sst1mpipe/resources/reference`` (``ref:<directory>``, made by the
scripts of ``scripts/reference``), with the quantities published only:

- ``crab_paper_2025/<mode>/ze<zenith>``: the SST-1M Crab paper (arXiv:2506.01733), ``tel1`` and
  ``strict_stereo`` (paper: stereo) at ``ze30``, ``ze50`` (energy bias and resolution, angular
  resolution, sensitivity), ``tel1``, ``tel2``, ``strict_stereo`` at ``ze20`` (effective area)
- sensitivity of other instruments: ``lhaaso_1yr``, ``hawc_1yr``, ``veritas_50h``, ``hess_50h``,
  ``astri_50h``

.. code-block:: console

    sst1mpipe-plot-irfs --irf strict_stereo=irf/strict_stereo \
        --reference "paper, 30 deg=ref:crab_paper_2025/strict_stereo/ze30" \
        --reference "ASTRI Mini-Array=ref:astri_50h" -o plots_paper

``sst1mpipe-time-to-detection`` gives the time to detect point sources at 5 sigma with the point-like
IRFs (``irf_point_like.fits.gz``) above minimum reconstructed energies: the spectrum of the source
(``sst1mpipe.performance.spectra``, ``--list-sources``) folded with the effective area and the energy
dispersion, the background of the IRFs in the ON region (``RAD_MAX``), 5 OFF regions, and the conditions
of the sensitivity of ctapipe (at least 10 signal events, signal at least 5 % of the background):

.. code-block:: console

    sst1mpipe-time-to-detection --irf irf/optimized/stereo/ze20 --source CRAB_HEGRA \
        --source MRK421_VERITAS_LOW --energy-min 1 --energy-min 3 -o time_to_detection.png
