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

``sst1mpipe-apply-models`` applies the models of the node of a production of models closest to the
pointing of the DL1 file (mono DL1 of ``sst1mpipe-process``, or stereo DL1 of
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

To make IRFs from MC DL2 files, one can run ``sst1mpipe_mc_make_irfs`` script, which currently produces only full enclosure IRFs, so it has to be provided with 
diffuse protons and diffuse gammas. The script applies event selection cuts defined in the config file (``event_selection``), including cut on gammaness. 
The gammaness cut can be either global (one number independent on energy) or energy dependent (gammaness distribution naturaly depends on energy, so using optimised 
energy dependent gammaness cut results in performance improvement). Global gammaness cut can be set in the config file (``global_gammaness_cut`` field), while energy dependent
cuts, must be provided as ``HDF5`` table, where the cuts for individual energy bins are stored, using parameter 
``--gammaness-cut-dir``. These tables can be generated with ``sst1mpipe_mc_performance`` (see :ref:`performance`). Pre-calculated energy dependent gammaness 
cuts are stored on Calculus for mono/stereo and different zenith angles:

``/data/work/analysis/MC/prod_january_2023/$SST1MPIPE_VER/performance/*_performance_*``

Again, one should make sure that the event selection applied to produce the cuts is the same as for IRFs. The IRF maker creates 
some directory structure inside the ``--output-dir``, automaticaly recognizing proper bin in zenith, azimuth, NSB level and gammaness cut applied. 
This directory structure should remain untouched for :ref:`dl2_dl3` to work properly.

Output IRF files are fully compatible with gammapy and may be read and explored with the use of gammapy funkcionalities:

.. code-block:: console

    from gammapy.irf import (
        EffectiveAreaTable2D,
        PSF3D,
        EnergyDispersion2D
        Background2D
    )
    irf_filename = 'SST1M_tel_021_Zen30deg_gcut0.75_irfs.fits'
    aeff = EffectiveAreaTable2D.read(irf_filename, hdu="EFFECTIVE AREA")
    edisp = EnergyDispersion2D.read(irf_filename, hdu="ENERGY DISPERSION")
    psf = PSF3D.read(irf_filename, hdu='POINT SPREAD FUNCTION')
    bg_2d = Background2D.read(irf_filename, hdu='BACKGROUND')


.. _dl2_dl3:

DL2 to DL3
~~~~~~~~~~

``sst1mpipe_data_dl2_dl3`` is a tool to create DL3 data files from DL2 data files. It is supposed to be provided with a directory with input DL2 files 
(typicaly a directory with DL2 for one source observed in one night, but can be run on larger sample as well). It merges the DL2 ``HDF5`` per-run 
files into per-wobble DL3 ``fits`` files containing only photon lists. It also finds proper IRF based on zenith, azimuth and NSB level for each input DL2 file. It creates 
per-night index files needed for further analysis in gammapy. The script applies event selection cuts defined in the config file (``event_selection``), 
including cut on gammaness, so **one should make sure that these are the same as used on MC to produce IRFs**. Energy dependent gammaness cuts can be used as well, following 
the same rules described in :ref:`dl2_irfs`.


.. _higher_level:

High level analysis
~~~~~~~~~~~~~~~~~~~

Output DL3 files produced with ``sst1mpipe_data_dl2_dl3`` are fully compatible with gammapy and may be further analyzed using gammapy tools. See e.g. 

* `Tutorial on 1D spectral analysis <https://docs.gammapy.org/1.2/tutorials/analysis-1d/spectral_analysis.html>`_

* `Tutorial on 2D ring background map <https://docs.gammapy.org/1.2/tutorials/analysis-2d/ring_background.html>`_

A typical use case is to run joint gammapy analysis on data from several nights. In such case one has to run ``create_hdu_indexes`` script to create 
HDU index files indexing all DL3s to be used in the final analysis.


.. _performance:

RF performance and sensitivity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TBD