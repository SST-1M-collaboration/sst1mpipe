.. _sst1m_analysis_workflow:

SST-1M analysis workflow
========================

.. _automatic_processing:

Automatic data processing
-------------------------

Raw data of both telescopes are processed up to DL3 every morning on **Calculus** by the scripts of
``scripts/daily`` (see ``scripts/daily/README.md``), run by cron (``scripts/daily/crontab``) on the
night before:

.. code-block:: console

    DATE=20260120 bash scripts/daily/run_daily.bash

The steps run on one machine, with all its cores. The same steps run with SLURM on camk with
``scripts/daily/submit_data.bash``. The products are in ``/data/work/analysis/Daily_analysis/daily_data/<DATE>/``
(``lists``, ``logs``, ``monitoring``, ``dl1``, ``dl2``, ``dl3``). Each data level has a directory per
telescope (``cs1``, ``cs2``) and one for the coincident events (``stereo``). The logs of each step and
raw file are in ``logs/``.

The models and IRFs are those of the simulations, processed on camk with ``scripts/camk/submit_mc.bash``
and copied to ``/data/work/analysis/MC/``. The settings are in ``scripts/daily/daily.env``.

DL3 products can be further analyzed using standard tools 
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

    .. include:: ../sst1mpipe/resources/config/sst1mpipe_mc_config_low_nsb.yaml
       :code: yaml

.. toggle:: 

    .. include:: ../sst1mpipe/resources/config/sst1mpipe_mc_config_high_nsb.yaml
       :code: yaml

Default config file for data:

.. toggle:: 

    .. include:: ../sst1mpipe/resources/config/sst1mpipe_data_config.yaml
       :code: yaml


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
  * ``pde_drop_factor`` - correction of the PDE drop due to the NSB in the simulations: the charges are divided by this factor (``factor`` of the camera calibration coefficients), which must match the PDE file of the simulation (``mc_pde_correction_factors.json``). ``null`` (real telescopes 21 and 22) applies no correction. ``sst1mpipe_mc_config_low_nsb.yaml`` and ``sst1mpipe_mc_config_high_nsb.yaml`` have the factors of the low and high NSB simulations

  The corrections which need the pedestal events are the camera calibration coefficients computed by ``sst1mpipe-monitoring`` (``CameraCoefficientsCalculator`` of the monitoring config, see :ref:`sst1mpipe_monitoring`): ``voltage_drop_correction`` (``none``, ``global`` or ``pixelwise`` correction of the voltage drop due to the NSB: the charges are divided by the charge drop of the ``SiPMVoltageDropModel``, see :ref:`nsb_level`) and ``flag_dead_pixels`` (the pixels without statistics of the ADC samples of the pedestal events are flagged as outliers and interpolated).

  The former ``telescope_calibration``, ``NsbCalibrator``, ``analysis``, ``event_selection``,
  ``telescope_equivalent_focal_length`` and ``window_transmittance`` settings are no longer in the configuration
  files (sst1mpipe < 0.9). The equivalent focal length is the one of the subarray description.


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

The IRFs are made with the tools of ctapipe, through their sst1mpipe versions, from the DL2 of
``sst1mpipe-apply-models`` (point gammas, diffuse gammas for the full enclosure IRFs, diffuse protons).
The sst1mpipe versions weight the events with the spectra of the SST-1M Crab paper (arXiv:2506.01733):
the Crab Nebula of MAGIC (``--crab MAGIC``, default; ``HEGRA`` as ctapipe), and the protons weighted to
protons + helium of DAMPE (``--background p_he``, default; ``protons`` as ctapipe). The relative
sensitivity (% Crab) is relative to the Crab of MAGIC.

1. ``sst1mpipe-optimize-event-selection``: gammaness and direction (theta) cuts, written to ``cuts.fits``:

   - ``sst1mpipe_optimize_cuts.yaml``: the cuts of the best sensitivity per reconstructed energy bin;
   - with ``sst1mpipe_irf_paper_cuts.yaml``: the cuts of the paper, a gammaness cut keeping 60 % of
     the gammas in each energy bin and one theta cut for all the energies keeping 80 % of the gammas
     (``GlobalThetaPercentileCuts``).
   - ``--EventSelectionOptimizer.optimization_algorithm=FixedThetaCuts``: the same gammaness cut and a
     fixed theta cut for all the energies, ``--FixedThetaCuts.theta_cut`` (0.2 deg: the cut of the analysis
     of the paper with sst1mpipe < 0.9).

   The stereo analysis has the events with both telescopes only: the stereo models predict only for them.
2. ``sst1mpipe-compute-irf`` (``sst1mpipe_compute_irf.yaml``): the IRFs with these cuts (effective area,
   energy dispersion, PSF or ``RAD_MAX``, background) in the GADF format of gammapy, and the
   benchmarks (energy bias and resolution, angular resolution, sensitivity in 50 h).
   Point-like IRFs: ``--spatial-selection-applied --DefaultFoVOffsetBins.fov_offset_n_bins=1``.
   The background of the sensitivity is the density of the gamma-like protons between
   ``Sensitivity2dMaker.fov_offset_min`` and ``Sensitivity2dMaker.fov_offset_max`` from the camera center,
   around the on-axis point source (as sst1mpipe < 0.9 and the SST-1M Crab paper), scaled to the size of the
   ON region; the IRFs keep the offset bins of ``DefaultFoVOffsetBins``. The sensitivity of an energy bin is
   computed only with at least ``MinimumStatisticsSensitivity2dMaker.min_simulated_signal`` simulated gammas
   and ``MinimumStatisticsSensitivity2dMaker.min_simulated_background`` simulated protons left after the cuts:

   .. literalinclude:: ../sst1mpipe/resources/config/sst1mpipe_compute_irf.yaml
      :language: yaml
      :start-at: Sensitivity2dMaker:
      :end-before: BackgroundRate2dMaker:

.. code-block:: console

    sst1mpipe-optimize-event-selection --config sst1mpipe_optimize_cuts.yaml \
        --config sst1mpipe_irf_paper_cuts.yaml -g gamma_point.dl2.h5 -p proton_test.dl2.h5 -o cuts.fits
    sst1mpipe-compute-irf --config sst1mpipe_compute_irf.yaml --cuts cuts.fits \
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

``sst1mpipe-dl2-dl3`` writes the DL3 (GADF format of gammapy) of a DL2 file of ``sst1mpipe-apply-models``,
with the cuts and the point-like IRFs of :ref:`dl2_irfs` (``--irf-directory``: a directory with
``cuts.fits`` and ``irf_point_like.fits.gz``, or with such directories ``ze<zenith>``: the zenith angle
closest to the mean pointing of each observation is used):

.. code-block:: console

    sst1mpipe-dl2-dl3 -i stereo/run.dl2.h5 --irf-directory irf/optimized/stereo -o dl3/

- one DL3 file per observation (``SST1M_<obs_id>.dl3.fits.gz``): ``EVENTS`` (the events passing the
  quality criteria and the gammaness cut of their reconstructed energy bin: reconstructed energy,
  RA/Dec, alt/az, gammaness), ``GTI``, ``POINTING``, and the IRFs (``EFFECTIVE AREA``,
  ``ENERGY DISPERSION``, ``RAD_MAX``, ``BACKGROUND``). The direction cut is not applied: the ``RAD_MAX``
  of the IRFs is the ON region of the point-like analysis of gammapy.
- the pointing of the telescopes (monitoring data), or of the target (warning) when the DL2 has no
  monitoring data
- the index files of gammapy (``hdu-index.fits.gz``, ``obs-index.fits.gz``) of all the DL3 files of the
  output directory (``--no-index`` when several runs are processed in parallel in the same directory,
  then ``create_hdu_indexes``)

The models, the cuts and the IRFs must be the same (the same mode: ``tel1`` / ``tel2`` for the mono
analysis of each telescope, ``stereo``). The dead time is not corrected yet (``DEADC = 1``), and the GTI
of an observation is from its first to its last event (issue #252).

The files follow the `Gamma Astro Data Formats <https://gamma-astro-data-formats.readthedocs.io>`_
(GADF, ``HDUCLASS = GADF``, ``HDUVERS = 0.3``), read by gammapy:

- ``EVENTS``: the mandatory columns ``EVENT_ID``, ``TIME`` (float64 s since ``MJDREFI``/``MJDREFF``,
  ``TIMESYS = TT``), ``RA``, ``DEC``, ``ENERGY``, and ``ALT``, ``AZ``, ``GAMMANESS``, ``MULTIP`` (telescopes
  of the direction); the mandatory keywords (``OBS_ID``, ``TSTART``, ``TSTOP``, ``ONTIME``, ``LIVETIME``,
  ``DEADC``, ``OBS_MODE = POINTING``, ``RA_PNT``, ``DEC_PNT``, ``EQUINOX``, ``RADECSYS``, ``ORIGIN``,
  ``TELESCOP = SST-1M``, ``INSTRUME`` = ``SST-1M-1``, ``SST-1M-2`` or ``SST-1M-1+2``, ``CREATOR``), the
  location of the array (``GEOLON``, ``GEOLAT``, ``ALTITUDE``) and ``OBJECT``, ``RA_OBJ``, ``DEC_OBJ``;
- ``GTI`` (``START``, ``STOP``), ``POINTING`` (``TIME``, ``RA_PNT``, ``DEC_PNT``, ``ALT_PNT``, ``AZ_PNT``);
- the point-like IRFs of ``ctapipe-compute-irf`` (``HDUCLAS3 = POINT-LIKE``: ``AEFF_2D``, ``EDISP_2D``,
  ``RAD_MAX_2D``; ``BKG_2D``);
- the index files: ``HDU_INDEX`` (the HDU types of GADF: ``events``, ``gti``, ``aeff``, ``edisp``,
  ``rad_max``, ``bkg``; the ``POINTING`` HDU is not indexed), ``OBS_INDEX`` (``OBS_ID``, ``OBS_MODE``,
  ``RA_PNT``, ``DEC_PNT``, ``TSTART``, ``TSTOP`` and ``ALT_PNT``, ``AZ_PNT``, ``ZEN_PNT``, ``ONTIME``,
  ``LIVETIME``, ``DEADC``, dates, ``OBJECT``, ``TMID``, ``EVENT_COUNT``, ...).

``sst1mpipe.io.dl3.check_dl3_file`` and ``check_index`` list what a DL3 file or the index files of a
directory miss of the specification:

.. code-block:: python

    from sst1mpipe.io.dl3 import check_dl3_file, check_index

    print(check_dl3_file("dl3/SST1M_20240110210012.dl3.fits.gz"))  # [] if it follows GADF
    print(check_index("dl3/"))


.. _higher_level:

High level analysis
~~~~~~~~~~~~~~~~~~~

The DL3 files of ``sst1mpipe-dl2-dl3`` are analysed with gammapy (point-like IRFs:
``required_irf="point-like"``):

.. code-block:: python

    from gammapy.data import DataStore

    observations = DataStore.from_dir("dl3/").get_observations(required_irf="point-like")

See e.g.

* `Tutorial on 1D spectral analysis <https://docs.gammapy.org/1.2/tutorials/analysis-1d/spectral_analysis.html>`_

* `Tutorial on 2D ring background map <https://docs.gammapy.org/1.2/tutorials/analysis-2d/ring_background.html>`_

A joint analysis of several nights needs the index files of all their DL3 files:

.. code-block:: console

    create_hdu_indexes --input-files '/dl3/20260120/*.dl3.fits.gz /dl3/20260121/*.dl3.fits.gz' -o crab/


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
  ``stereo`` at ``ze30``, ``ze50`` (energy bias and resolution, angular
  resolution, sensitivity), ``tel1``, ``tel2``, ``stereo`` at ``ze20`` (effective area)
- sensitivity of other instruments: ``lhaaso_1yr``, ``hawc_1yr``, ``veritas_50h``, ``hess_50h``,
  ``astri_50h``

.. code-block:: console

    sst1mpipe-plot-irfs --irf stereo=irf/stereo \
        --reference "paper, 30 deg=ref:crab_paper_2025/stereo/ze30" \
        --reference "ASTRI Mini-Array=ref:astri_50h" -o plots_paper
