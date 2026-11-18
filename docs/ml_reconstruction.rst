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
``sst1mpipe-apply-models``               ``ctapipe-apply-models``                   models of the closest node of a production of models
=======================================  =========================================  ===================================================

The stereo direction of the MARS analysis (``StereoMARSCombiner``) is available in all of them.

Input: DL1 with the telescope frame parameters
----------------------------------------------

The disp reconstructor and the intensity weights of the stereo combination use the Hillas
parameters in the telescope frame (``hillas_fov_lon``, ``hillas_fov_lat``, ``hillas_psi``,
``hillas_intensity``). The MC configurations of sst1mpipe write the camera frame parameters:

.. code-block:: console

    sst1mpipe-process -i gamma.simtel.gz -o gamma.dl1.h5 --config sst1mpipe_mc_config_low_nsb.json \
        --ImageProcessor.use_telescope_frame=True

The training files are merged with ``ctapipe-merge``.

Training
--------

The configurations ``sst1mpipe/resources/config/sst1mpipe_train_<model>.yaml``
(``sst1mpipe.resources.TRAIN_CONFIG_FILES``) have the random forests, features and event selection
of ``sst1mpipe_mc_config_low_nsb.json``:

.. code-block:: console

    sst1mpipe-train-energy-regressor -i gamma_diffuse_train.dl1.h5 -o models/ze20_az0_nsb100/energy_regressor.pkl \
        --config sst1mpipe_train_energy_regressor.yaml --n-events 700000
    sst1mpipe-train-particle-classifier --signal gamma_diffuse_train.dl1.h5 --background proton_train.dl1.h5 \
        -o models/ze20_az0_nsb100/particle_classifier.pkl --config sst1mpipe_train_particle_classifier.yaml \
        --n-events 700000 --signal-fraction 0.5
    sst1mpipe-train-disp-reconstructor -i gamma_diffuse_train.dl1.h5 -o models/ze20_az0_nsb100/disp_reconstructor.pkl \
        --config sst1mpipe_train_disp_reconstructor.yaml --n-events 700000

One model is trained per telescope type (the two SST-1M telescopes are two types in the subarray).
The models are stereo models if the stereo features are in the configuration
(``HillasReconstructor_tel_impact_distance``, ``HillasReconstructor_h_max``, commented out in the
configurations).

After the training, the plots of the models are made in ``--plot-dir`` (by default
``<model name>_plots`` next to the model, e.g. ``models/ze20_az0_nsb100/energy_regressor_plots``),
with ctaplot:

- ``<model name>_feature_importance.png``: feature importances of the model of each telescope type
  (the norm regressor and the sign classifier of the disp);
- ``<model name>_performance.png``: performance from the cross validation of the training tools
  (``CrossValidator.n_cross_validations``, 5 by default; its predictions are written in
  ``<model name>_cv.h5`` in the plot directory if ``--cv-output`` is not given).

The cross validation tables are per telescope type name: the two SST-1M telescopes, of the same
type name, are in the same performance plot.

Application: the closest models
-------------------------------

A production of models has one subdirectory per node ``zeXX_azXX_nsbXX`` (zenith and azimuth in
degrees, NSB rate in MHz) with the models of the node (``*.pkl``). ``sst1mpipe-apply-models``
applies the models of the node closest to the pointing of the input file: the closest NSB rate if
``--nsb-rate`` is given, then the closest zenith, then the closest azimuth (median pointing of the
telescope events):

.. code-block:: console

    sst1mpipe-apply-models -i run.dl1.h5 --models-directory models -o run.dl2.h5

With ``--reconstructor`` (``-r``), the models given are used instead.

Stereo direction of the MARS analysis
-------------------------------------

``StereoMARSCombiner`` is selected in the configuration of the disp reconstructor (as in
``sst1mpipe_train_disp_reconstructor.yaml``):

.. code-block:: yaml

    DispReconstructor:
      stereo_combiner_cls: StereoMARSCombiner
      StereoMARSCombiner:
        weights: intensity

Each telescope gives two possible source positions: at the reconstructed disp norm on either side
of the image centroid along its major axis. The combination of signs whose positions are the
closest to each other (smallest sum of the pairwise angular distances; 4 combinations for 2
telescopes) is chosen, and the stereo direction is the weighted mean of its positions. The disp sign
predicted by the classifier is only used for the events of a single telescope. The
``disp_goodness_of_fit`` of the stereo events is the sum of the angular distances of the chosen
positions (deg). The energy and the gammaness are combined as by the ``StereoMeanCombiner``.

The combiner is a class of sst1mpipe (``sst1mpipe.ml``): the models using it are applied with the
sst1mpipe tools (or any program importing ``sst1mpipe.ml``), not with ``ctapipe-apply-models``.

Differences with the random forests of sst1mpipe
------------------------------------------------

==============================================  =================================================================
same                                            random forests and their settings, features, event selection,
                                                number of training events, gamma/proton fraction, one model per
                                                telescope, intensity weighted stereo mean of the energy and
                                                gammaness, MARS stereo direction
different                                       energy trained on ln(E) (log10(E) in sst1mpipe), one feature list
                                                for the disp norm and sign, DL2 columns and tables of ctapipe, one
                                                file per model type (with its configuration)
not (yet) with the ctapipe reconstructors       ``disp_vector``, misdirection model and event classes, reconstructed
                                                energy and disp as gammaness features (by hand: apply the energy and
                                                disp models to a separate classifier training sample)
new                                             cross validation, training by chunks (larger samples)
==============================================  =================================================================
