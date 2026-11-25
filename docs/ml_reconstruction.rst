DL1 to DL2 with the ctapipe reconstructors
==========================================

The energy, the gamma/hadron separation (gammaness) and the direction of the events can be
reconstructed with the machine learning reconstructors of ctapipe (``ctapipe.reco``), trained and
applied by sst1mpipe tools built on the ctapipe tools. They have the options and configuration
sections of the ctapipe tools, and write the DL2 of ctapipe (tables and column names, e.g.
``RandomForestRegressor_energy``, ``RandomForestClassifier_prediction``, ``disp_alt``):

=======================================  =========================================  ===================================================
sst1mpipe tool                           ctapipe tool                               added by sst1mpipe
=======================================  =========================================  ===================================================
``sst1mpipe-train-energy-regressor``     ``ctapipe-train-energy-regressor``         plots: feature importances, energy resolution, bias
                                                                                    and migration matrix
``sst1mpipe-train-particle-classifier``  ``ctapipe-train-particle-classifier``      plots: feature importances, gammaness, ROC curves
                                                                                    (all energies and per true energy)
``sst1mpipe-train-disp-reconstructor``   ``ctapipe-train-disp-reconstructor``       plots: feature importances, disp norm, sign accuracy
``sst1mpipe-apply-models``               ``ctapipe-apply-models``                   pointing of the target if no telescope pointing,
                                                                                    stereo combiner of the direction
=======================================  =========================================  ===================================================

Input: DL1 with the telescope frame parameters
----------------------------------------------

The disp reconstructor and the intensity weights of the stereo combination use the Hillas
parameters in the telescope frame (``hillas_fov_lon``, ``hillas_fov_lat``, ``hillas_psi``,
``hillas_intensity``), written by the sst1mpipe configurations (``ImageProcessor.use_telescope_frame``):

.. code-block:: console

    sst1mpipe-process -i gamma.simtel.gz -o gamma.dl1.h5 --config sst1mpipe_mc_config_low_nsb.json

The training files are merged with ``ctapipe-merge``.

Training
--------

The configurations ``sst1mpipe/resources/config/sst1mpipe_train_<model>.yaml``
(``sst1mpipe.resources.TRAIN_CONFIG_FILES``) have the random forests and features of the former
``sst1mpipe_mc_config_low_nsb.json``, and the event selection of the SST-1M Crab paper
(arXiv:2506.01733): intensity > 45 p.e., ``leakage_intensity_width_2`` < 0.7 (``MLQualityQuery``,
applied at the training and at the application of the models: the other telescope events have no
DL2):

.. code-block:: console

    sst1mpipe-train-energy-regressor -i gamma_diffuse_train.dl1.h5 -o models/nsb100/energy_regressor.pkl \
        --config sst1mpipe_train_energy_regressor.yaml --n-events 700000
    sst1mpipe-train-particle-classifier --signal gamma_diffuse_train.dl1.h5 --background proton_train.dl1.h5 \
        -o models/nsb100/particle_classifier.pkl --config sst1mpipe_train_particle_classifier.yaml \
        --n-events 700000 --signal-fraction 0.5
    sst1mpipe-train-disp-reconstructor -i gamma_diffuse_train.dl1.h5 -o models/nsb100/disp_reconstructor.pkl \
        --config sst1mpipe_train_disp_reconstructor.yaml --n-events 700000

ctapipe trains one model per telescope type: the two SST-1M telescopes are two types in the subarray
(same name ``SST_1M_DigiCam``), so one model per telescope.
Two sets of models, as in the SST-1M Crab paper:

- **mono models** (``sst1mpipe_train_<model>.yaml``): the image parameters only, for the analysis of
  one telescope;
- **stereo models** (``sst1mpipe_train_<model>_stereo.yaml``, ``sst1mpipe.resources.STEREO_TRAIN_CONFIG_FILES``):
  in addition the stereo parameters of the HillasReconstructor, ``HillasReconstructor_tel_impact_distance``
  (impact distance of the telescope) and ``HillasReconstructor_h_max`` (height of the shower maximum),
  for the energy, the gammaness, and the disp norm and sign. They are trained on the events reconstructed
  by both telescopes (quality criterion ``stereo``) and give predictions for these events only. Still one
  model per telescope, whose predictions are combined (intensity weighted mean of the energy and the
  gammaness, stereo direction of the paper: ``StereoDISPSignCombiner``). The simulations have the stereo parameters (the HillasReconstructor runs in
  ``sst1mpipe-process``), the stereo data those of ``sst1mpipe-stereo``.

After the training, the plots of the models are made in ``--plot-dir`` (by default
``<model name>_plots`` next to the model, e.g. ``models/nsb100/energy_regressor_plots``),
with ctaplot:

- one figure per quantity and per telescope, ``<model name>_tel<id>_<quantity>.png``: feature
  importances, and the performance from the cross validation of ctapipe
  (``CrossValidator.n_cross_validations``, 5 by default; written in ``<model name>_cv.h5`` in the
  plot directory if ``--cv-output`` is not given): energy resolution, bias and migration; gammaness
  and ROC curves; disp norm migration, resolution, bias and sign accuracy;
- the same axis limits for all the figures (``sst1mpipe.ml.plots.LIMITS``: relative errors within
  0.5), bias and resolution as the benchmarks of ctapipe (median and half of the 16-84 % interval of
  the relative error, energy bins with at least 10 events).

The stereo performance is the one of the DL2 (benchmarks of ``ctapipe-compute-irf``).

Some plots of the models trained on the three test simtel files (``tests``, a few hundred events:
the statistics are low), telescope 1:

.. image:: _static/ml/energy_regressor_tel1_feature_importance.png
   :width: 49%
.. image:: _static/ml/energy_regressor_tel1_resolution.png
   :width: 49%
.. image:: _static/ml/energy_regressor_tel1_bias.png
   :width: 49%
.. image:: _static/ml/energy_regressor_tel1_migration.png
   :width: 49%
.. image:: _static/ml/particle_classifier_tel1_gammaness.png
   :width: 49%
.. image:: _static/ml/particle_classifier_tel1_roc.png
   :width: 49%
.. image:: _static/ml/disp_reconstructor_tel1_disp_resolution.png
   :width: 49%
.. image:: _static/ml/disp_reconstructor_tel1_disp_sign.png
   :width: 49%

Zenith angle and NSB
--------------------

The models of an NSB level are trained on the training files of all its simulated zenith angles
merged, with the pointing of the telescope as feature (``cos_zenith``, from
``telescope_pointing_altitude``; not the azimuth, the reconstruction is not sensitive to it): the
same models for all the pointings, in ``models/nsb<rate>/``. The NSB level is chosen with the
models given to ``sst1mpipe-apply-models``.

A random forest does not extrapolate: the events of a zenith angle between two simulated ones are
reconstructed close to the closest simulated one, and outside the simulated zenith angles as the
closest one (e.g. 40 deg as 30 deg with simulations at 10, 20, 30 deg). The simulated zenith
angles must cover the observations.

Application
-----------

.. code-block:: console

    sst1mpipe-apply-models -i run.dl1.h5 -o run.dl2.h5 -r models/nsb100/energy_regressor.pkl \
        -r models/nsb100/particle_classifier.pkl -r models/nsb100/disp_reconstructor.pkl

The pointing of the telescopes is in the simulations, and in the data processed with the monitoring
data (interpolated at the time of each event). Without it (data processed without monitoring
data), the pointing of the target (RA/Dec of the run) transformed into alt/az at the time of each
event is used, with a warning.

Stereo direction of the SST-1M Crab paper
-----------------------------------------

The stereo disp models use ``StereoDISPSignCombiner`` (``sst1mpipe_train_disp_reconstructor_stereo.yaml``),
the stereo direction of the SST-1M Crab paper (arXiv:2506.01733, Sec. 4.3.2, the approach of MAGIC, as
``stereo_reconstruction`` of sst1mpipe < 0.9):

.. code-block:: yaml

    DispReconstructor:
      stereo_combiner_cls: StereoDISPSignCombiner
      StereoDISPSignCombiner:
        weights: intensity
        min_telescopes: 2

- only the disp norm of each telescope is used, not its sign classifier: each telescope gives two
  possible source positions, at its disp norm on either side of the image centroid along the major axis;
- of the combinations of the signs of the telescopes (4 for two telescopes), the one whose positions are
  the closest is chosen (``disp_sign_combination``: smallest angular distance, for more telescopes the
  smallest sum of the distances of the pairs);
- the direction is the weighted mean of these positions (``weights``, as ``StereoMeanCombiner``:
  intensity, as the paper);
- ``sign_method: axes_crossing``: the former sst1mpipe method instead (each telescope takes the sign
  whose position is the closest to the crossing of the major axes, ``axes_crossing_disp_signs``), for
  comparison;
- stereo events are the events with at least ``min_telescopes`` valid telescopes, i.e. passing the
  quality cuts of the models (paper: both telescopes, intensity > 45 p.e., leakage2 < 0.7). The other
  ones get an invalid direction, so they are not in the IRFs. ``min_telescopes: 1`` keeps the events
  with one telescope, with the sign of its classifier.

The energy and the gammaness are combined by the ``StereoMeanCombiner`` (intensity weighted mean). The
combiner is a class of sst1mpipe (``sst1mpipe.ml``): the models using it are applied with the sst1mpipe
tools (or any program importing ``sst1mpipe.ml``), not with ``ctapipe-apply-models``.

The combiner is saved with the disp model at its training. ``sst1mpipe-apply-models`` can apply the
same models with another combiner, e.g. to compare the stereo directions:

.. code-block:: console

    sst1mpipe-apply-models -i gamma_point.dl1.h5 -o gamma_point.dl2.h5 -r models/stereo/disp_reconstructor.pkl ... \
        --ApplyModels.disp_stereo_combiner=StereoDISPSignCombiner --StereoDISPSignCombiner.min_telescopes=2

Differences with the random forests of sst1mpipe
------------------------------------------------

==============================================  =================================================================
same                                            random forests and their settings, features, event selection,
                                                number of training events, gamma/proton fraction, one model per
                                                telescope, intensity weighted stereo mean of the energy and
                                                gammaness, stereo direction (``StereoDISPSignCombiner``)
different                                       energy trained on ln(E) (log10(E) in sst1mpipe), one feature list
                                                for the disp norm and sign, DL2 columns and tables of ctapipe, one
                                                file per model type (with its configuration)
not (yet) with the ctapipe reconstructors       ``disp_vector``, misdirection model and event classes, reconstructed
                                                energy and disp as gammaness features (by hand: apply the energy and
                                                disp models to a separate classifier training sample)
new                                             cross validation, training by chunks (larger samples)
==============================================  =================================================================

API
---

.. autoclass:: sst1mpipe.ml.StereoDISPSignCombiner

.. autofunction:: sst1mpipe.ml.disp_sign_combination

.. autofunction:: sst1mpipe.ml.axes_crossing_disp_signs
