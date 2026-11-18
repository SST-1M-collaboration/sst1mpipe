import logging
import os

import astropy.units as u
import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import AltAz
from astropy.table import QTable
from astropy.time import Time
from ctaplot.ana import angular_separation_altaz, logbin_mean
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import train_test_split

from sst1mpipe.analysis import add_reco_ra_dec
from sst1mpipe.io import (
    check_outdir,
)
from sst1mpipe.utils import (
    camera_to_altaz,
    disp_to_pos,
    disp_vector,
    get_event_sample,
    get_finite,
    SIMULATION_OBSTIME,
    telescope_location,
    get_tel_id,
    mix_gamma_proton,
    remove_stereo,
)


def plot_feature_importance(
        model, features=None,
        outfile=None,
        telescope=None):
    """
    Save figures with feature importances for given RF model

    Parameters
    ----------
    model: sklearn.ensemble.RandomForestClassifier or sklearn.ensemble.RandomForestRegressor
    features: list of string
    outfile: str
        output file path
    telescope: str (tel_001/tel_002)

    Returns
    -------

    """

    plt.figure()
    ax = plt.gca()

    importances = model.feature_importances_
    std = np.std([tree.feature_importances_ for tree in model.estimators_], axis=0)
    indices = np.argsort(importances)

    ordered_features = [features[index] for index in indices]

    ax.set_title("Feature importances (impurity based), " + telescope)

    ax.barh(ordered_features,
            importances[indices],
            xerr=std[indices],
            align="center",
            )

    ax.set_yticks(range(len(features)))
    ax.set_yticklabels(ordered_features)
    ax.grid()
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)


def reco_source_position_sky(
        cog_x, cog_y, disp_dx, disp_dy,
        focal_length, pointing_alt, pointing_az,
        location, times=None):

    src_x, src_y = disp_to_pos(disp_dx, disp_dy, cog_x, cog_y)
    return camera_to_altaz(src_x, src_y, focal_length, pointing_alt, pointing_az, location=location, times=times)


def train_rf_energy(
        data, config=None, plot=False,
        outdir=None, telescope=None, stereo=False):
    """
    Train RF Regressor for energy reconstruction

    Parameters
    ----------
    data: pandas.DataFrame
    config: dict
    plot: bool
    outdir: str
        output file path
    telescope: str (tel_001/tel_002)
    stereo: bool

    Returns
    -------
    Trained model, sklearn.ensemble.RandomForestRegressor

    """

    logging.info("Training Random Forest Regressor for Energy Reconstruction...")
    energy_regression_args = config['random_forest_regressor_args']
    features = config['energy_regression_features']

    if not stereo:
        features = remove_stereo(features)

    logging.info("Given features: %s", features)
    logging.info("Number of events for training: %d", len(data))

    reg = RandomForestRegressor(**energy_regression_args)
    reg.fit(data[features],
            data['log_true_energy'])

    logging.info("Energy reconstruction model trained!")

    if plot:
        plot_feature_importance(reg, features=features, outfile=outdir + '/reg_energy_features_' + str(telescope) + '.png', telescope=telescope)

    joblib.dump(reg, outdir + '/reg_energy_' + str(telescope) + '.sav', compress=3)
    del reg


def train_disp_vector(
        data, config=None, plot=False,
        outdir=None, telescope=None, stereo=False):
    """
    Train RF Regressor for DISP reconstruction (vector of disp_dx, disp_dy)

    Parameters
    ----------
    data: pandas.DataFrame
    config: dict
    plot: bool
    outdir: str
        output file path
    telescope: str (tel_001/tel_002)
    stereo: bool

    Returns
    -------
    Trained model, sklearn.ensemble.RandomForestRegressor

    """

    logging.info("Training Random Forest Regressor for DISP Reconstruction...")
    predict_features = ['disp_dx', 'disp_dy']
    disp_regression_args = config['random_forest_regressor_args']
    features = config['disp_regression_features']

    if not stereo:
        features = remove_stereo(features)

    logging.info("Given features: %s", features)
    logging.info("Number of events for training: %d", len(data))

    reg = RandomForestRegressor(**disp_regression_args)
    x = data[features]
    y = np.transpose([data[f] for f in predict_features])
    reg.fit(x, y)

    logging.info("DISP reconstruction model trained!")

    if plot:
        plot_feature_importance(reg, features=features, outfile=outdir + '/reg_disp_features_' + str(telescope) + '.png', telescope=telescope)

    joblib.dump(reg, outdir + '/reg_disp_' + str(telescope) + '.sav', compress=3)
    del reg

def train_disp_norm(
        data, config=None, plot=False,
        outdir=None, telescope=None, stereo=False):
    """
    Train RF Regressor for DISP norm reconstruction

    Parameters
    ----------
    data: pandas.DataFrame
    config: dict
    plot: bool
    outdir: str
        output file path
    telescope: str (tel_001/tel_002)
    stereo: bool

    Returns
    -------
    Trained model, sklearn.ensemble.RandomForestRegressor

    """

    logging.info("Training Random Forest Regressor for DISP Reconstruction...")
    disp_regression_args = config['random_forest_regressor_args']
    features = config['disp_regression_features']

    if not stereo:
        features = remove_stereo(features)

    logging.info("Given features: %s", features)
    logging.info("Number of events for training: %d", len(data))

    reg = RandomForestRegressor(**disp_regression_args)
    reg.fit(data[features],
            data['disp_norm'])

    logging.info("DISP NORM reconstruction model trained!")

    if plot:
        plot_feature_importance(reg, features=features, outfile=outdir + '/reg_disp_norm_features_' + str(telescope) + '.png', telescope=telescope)

    joblib.dump(reg, outdir + '/reg_disp_norm_' + str(telescope) + '.sav', compress=3)
    del reg


def train_disp_sign(
        data, config=None, plot=False,
        outdir=None, telescope=None, stereo=False):
    """
    Train RF Classifier for DISP sign reconstruction

    Parameters
    ----------
    data: pandas.DataFrame
    config: dict
    plot: bool
    outdir: str
        output file path
    telescope: str (tel_001/tel_002)
    stereo: bool

    Returns
    -------
    Trained model, sklearn.ensemble.RandomForestClassifier

    """

    logging.info("Training Random Forest Classifier for DISP sign...")
    classification_args = config['random_forest_classifier_args']
    features = config["disp_classification_features"]

    if not stereo:
        features = remove_stereo(features)

    logging.info("Given features: %s", features)
    logging.info("Number of events for training: %d", len(data))

    clf = RandomForestClassifier(**classification_args)

    clf.fit(data[features],
            data['disp_sign'])
    logging.info("DISP SIGN classifier trained!")

    if plot:
        plot_feature_importance(clf, features=features, outfile=outdir + '/class_disp_sign_features_' + str(telescope) + '.png', telescope=telescope)

    joblib.dump(clf, outdir + '/class_disp_sign_' + str(telescope) + '.sav', compress=3)
    del clf


def train_rf_separation(
        data, config=None, plot=False,
        outdir=None, telescope=None, stereo=False):
    """
    Train RF Classifier for gamma/hadron separation

    Parameters
    ----------
    data: pandas.DataFrame
    config: dict
    plot: bool
    outdir: str
        output file path
    telescope: str (tel_001/tel_002)
    stereo: bool

    Returns
    -------
    Trained model, sklearn.ensemble.RandomForestClassifier

    """

    logging.info("Training Random Forest Classifier for Gamma/Hadron separation...")
    classification_args = config['random_forest_classifier_args']
    features = config["particle_classification_features"]

    if not stereo:
        features = remove_stereo(features)

    logging.info("Given features: %s", features)
    logging.info("Number of events for training: %d", len(data))

    clf = RandomForestClassifier(**classification_args)

    clf.fit(data[features],
            data['true_shower_primary_id'])  # 0 .. gamma, 101 .. proton
    logging.info("Gamma/hadron classifier trained!")

    if plot:
        plot_feature_importance(clf, features=features, outfile=outdir + '/class_gh_features_' + str(telescope) + '.png', telescope=telescope)

    joblib.dump(clf, outdir + '/class_gh_' + str(telescope) + '.sav', compress=3)
    del clf


def train_rf_misdirection(
    data, config=None,
    outdir=None, telescope=None,
    plot=None):

    regression_args = config['random_forest_regressor_args']
    if telescope == 'stereo':
        features = config['misdirection_regression_features_stereo']
    else:
        features = config['misdirection_regression_features_mono']

    reg = RandomForestRegressor(**regression_args)
    reg.fit(data[features], data['log_misdirection'])

    if plot:
        plot_feature_importance(reg, features=features, outfile=outdir + '/reg_mis_features_' + str(telescope) + '.png', telescope=telescope)

    joblib.dump(reg, outdir + '/reg_misdirection_' + str(telescope) + '.sav', compress=3)
    del reg


def train_models(
        params_gamma, params_protons, config=None,
        plot=None, outdir=None, telescope=None,
        stereo=False):
    """
    Run training of all RF models on MC gammas and protons

    Parameters
    ----------
    params_gamma: pandas.DataFrame
    params_protons: pandas.DataFrame
    config: dict
    plot: bool
    outdir: str
        output file path
    telescope: str (tel_001/tel_002)
    stereo: bool

    Returns
    -------

    """

    N_training_events = config["analysis"]["n_training_events"]
    gp_ratio = config["analysis"]["gamma_to_proton_training_ratio"]
    gamma_training_sample = get_event_sample(params_gamma, max_events=N_training_events)

    # Training models
    if (config['disp_method'] == 'disp_vector') & (not stereo):
        train_disp_vector(gamma_training_sample, config=config, plot=plot, outdir=outdir, telescope=telescope, stereo=stereo)
    elif (config['disp_method'] == 'disp_norm_sign') & (not stereo):
        train_disp_norm(gamma_training_sample, config=config, plot=plot, outdir=outdir, telescope=telescope, stereo=stereo)
        train_disp_sign(gamma_training_sample, config=config, plot=plot, outdir=outdir, telescope=telescope, stereo=stereo)
    elif stereo:
        train_disp_norm(gamma_training_sample, config=config, plot=plot, outdir=outdir, telescope=telescope, stereo=stereo)
        logging.info("STEREO reconstruction: ONLY DISP NORM REGRESSOR WILL BE TRAINED (NOT THE SIGN CLASSIFIER), regardless the method set in the config file.")
    else:
        raise Exception("ERROR: disp_method in the config file not recognized, use one of these: disp_vector, disp_norm_sign")

    train_rf_energy(gamma_training_sample, config=config, plot=plot, outdir=outdir, telescope=telescope, stereo=stereo)

    if params_protons is not None:

        # We want to use reco energy and disp_norm as features for g/h classifier (we skip disp_sign as it
        # would be quite painful to implement again for stereo). To do this, we need to split training gammas into two subsets,
        # train energy and disp reconstructors on one of them, reco these two parameters in second of them, and then use the
        # second one as training sample for g/h separator
        if 'log_reco_energy' in config['particle_classification_features']:

            #N_gammas_train = int(len(gamma_training_sample)/2.)
            gamma_train_temp, gamma_test_temp = train_test_split(gamma_training_sample, test_size=0.33, random_state=42)

            # Training temporary models
            outdir_tmp = outdir + '/tmp/'
            check_outdir(outdir_tmp)
            logging.info('Training temporary models for energy and disp_norm reco for training of g/h separation.')
            train_disp_norm(gamma_train_temp, config=config, plot=plot, outdir=outdir_tmp, telescope=telescope, stereo=stereo)
            train_rf_energy(gamma_train_temp, config=config, plot=plot, outdir=outdir_tmp, telescope=telescope, stereo=stereo)

            # Reco energy and disp
            energy_regression_features = config["energy_regression_features"]
            disp_regression_features = config["disp_regression_features"]

            if not stereo:
                energy_regression_features = remove_stereo(energy_regression_features)
                disp_regression_features = remove_stereo(disp_regression_features)

            energy_re = joblib.load(os.path.join(outdir_tmp, 'reg_energy_' + str(telescope) + '.sav'))
            logging.info('Loading model: ' + outdir_tmp + '/reg_energy_' + str(telescope) + '.sav')
            logging.info('Temporary energy reconstruction for g/h separation training..')
            gamma_test_temp['log_reco_energy'] = energy_re.predict(gamma_test_temp[energy_regression_features])
            params_protons['log_reco_energy'] = energy_re.predict(params_protons[energy_regression_features])

            disp_re = joblib.load(os.path.join(outdir_tmp, 'reg_disp_norm_' + str(telescope) + '.sav'))
            logging.info('Loading model: ' + outdir_tmp + '/reg_disp_norm_' + str(telescope) + '.sav')
            logging.info('Temporary DISP norm reconstruction for g/h separation training..')
            gamma_test_temp['reco_disp_norm'] = disp_re.predict(gamma_test_temp[disp_regression_features])
            params_protons['reco_disp_norm'] = disp_re.predict(params_protons[disp_regression_features])

            gh_training_dataset = mix_gamma_proton(gamma_test_temp, params_protons, max_events=N_training_events, gp_training_ratio=gp_ratio)
            gh_training_dataset = get_finite(gh_training_dataset, config=config, stereo=stereo)

            train_rf_separation(gh_training_dataset, config=config, plot=plot, outdir=outdir, telescope=telescope, stereo=stereo)
            logging.info('Gamma/hadron separator trained using reco_energy and disp_norm as additional features.')

        else:
            gh_training_dataset = mix_gamma_proton(params_gamma, params_protons, max_events=N_training_events, gp_training_ratio=gp_ratio)
            train_rf_separation(gh_training_dataset, config=config, plot=plot, outdir=outdir, telescope=telescope, stereo=stereo)
    else:
        logging.warning('RF classifier cannot be trained, dus to missing input DL1 proton file!')
    logging.info('All Random Forest models trained and stored here: %s', outdir)


def reco_misdirection(dl2, models_dir=None, config=None, telescope=None):

    if telescope == 'stereo':

        # In data, we need to rename columns _tel21 -> tel1, tel22 -> tel2
        for tt,ttn in zip(['_tel21', '_tel22'], ['_tel1', '_tel2'], strict=False):
            cols_to_rename = dl2.columns[dl2.columns.str.contains(tt)]
            rename_dict = {col: col.replace(tt, ttn) for col in cols_to_rename}
            dl2 = dl2.rename(columns=rename_dict)

        dl2['log_camera_frame_hillas_intensity_tel1'] = np.log10(dl2['camera_frame_hillas_intensity_tel1'])
        dl2['log_camera_frame_hillas_intensity_tel2'] = np.log10(dl2['camera_frame_hillas_intensity_tel2'])

    if telescope == 'stereo':
        features = config['misdirection_regression_features_stereo']
    else:
        features = config['misdirection_regression_features_mono']

    # Check finite
    mask = np.ones(len(dl2), dtype=bool)
    for key in features:
        try:
            mask &= np.isfinite(dl2[key])
        except KeyError:
            logging.warning(f'{key} column not in data.')
    dl2_finite = dl2[mask].copy()

    mis_mono_re = joblib.load(models_dir + '/reg_misdirection_'+telescope+'.sav')
    logging.info('Misdirection reconstruction..')
    dl2_finite['log_reco_misdirection'] = mis_mono_re.predict(dl2_finite[features])
    return dl2_finite


def get_evttype_edges(dl2, config=None, percentiles=None):
    if percentiles is None:
        percentiles = [25, 50, 75]

    ebins = np.logspace(config['analysis']["log_energy_min_tev"], config['analysis']["log_energy_max_tev"], config['analysis']["n_energy_bins"])
    energy_center = logbin_mean(ebins)

    edges = []
    for i in range(len(ebins)-1):
        mask = (dl2['reco_energy'] > ebins[i]) & (dl2['reco_energy'] < ebins[i+1])

        if sum(mask):
            evttype_edges = np.percentile(
                    dl2['log_reco_misdirection'][mask],
                    percentiles
                )
            evttype_edges = np.concatenate(
                ([-np.inf], evttype_edges, [np.inf])
            )
        else:
            evttype_edges = np.array([np.nan, np.nan, np.nan, np.nan, np.nan])
        edges.append(evttype_edges)
    edges = np.array(edges)


    cut_table = QTable()
    cut_table["low"] = ebins[:-1]
    cut_table["high"] = ebins[1:]
    cut_table["center"] = energy_center
    for i in range(len(percentiles)+2):
        key = 'log_midirection_edge_' + str(i)
        cut_table[key] = np.asanyarray(np.nan, edges.dtype)
        cut_table[key].value[:] = edges[:,i]

    return cut_table


def classify_evt_types(dl2, edges=None):

    dl2['event_type'] = 0
    ebins_low = edges["low"]
    ebins_high = edges["high"]

    evttype_edges = [s for s in edges.keys() if 'log_midirection_' in s]

    for i in range(len(ebins_low)):
        mask_e = (dl2['reco_energy'] > ebins_low[i]) & (dl2['reco_energy'] < ebins_high[i])

        for j in range(len(evttype_edges)-1):
            mask_class = (dl2['log_reco_misdirection'] > edges[evttype_edges[j]][i]) & (dl2['log_reco_misdirection'] <= edges[evttype_edges[j+1]][i])
            dl2.loc[mask_e & mask_class, 'event_type'] = j+1

    return dl2


def apply_models(
        dl1, models_dir=None, config=None,
        telescope=None, stereo=False, mc=True, subarray=None):
    """
    Applies all RF models from models directory to reconstruct
    all events from DL1 to DL2. Works on both MC/data, mono/stereo.
    For stereo reconstruction, only DISP norm is reconstructed,
    and both possible signs are stored in the output DataFrame for
    each event (MARS-like stereo reconstruction). All quantities
    are reconstructed per-telescope even in stereo. Final averaging
    for stereo reconstruction is done in .reco.stereo_reconstruction()

    Parameters
    ----------
    dl1: pandas.DataFrame
    models_dir: str
        models dir path
    config: dict
    telescope: str (tel_001/tel_002)
    stereo: bool
    mc: bool
    subarray: ctapipe.instrument.SubarrayDescription
        subarray of the input file (``SubarrayDescription.from_hdf(input_file)``), for the telescope position

    Returns
    -------
    pandas.DataFrame

    """
    if subarray is None:
        raise ValueError("subarray is required (SubarrayDescription.from_hdf(input_file)): the telescope positions")

    if not mc:
        if '1' in telescope:
            telescope_model = 'tel_001'
        else:
            telescope_model = 'tel_002'
    else:
        telescope_model = telescope

    dl2 = dl1.copy()

    energy_regression_features = config["energy_regression_features"]
    gh_clas_features = config["particle_classification_features"]
    disp_regression_features = config["disp_regression_features"]
    disp_classification_features = config["disp_classification_features"]

    if not stereo:
        energy_regression_features = remove_stereo(energy_regression_features)
        gh_clas_features = remove_stereo(gh_clas_features)
        disp_regression_features = remove_stereo(disp_regression_features)
        disp_classification_features = remove_stereo(disp_classification_features)

    try:
        energy_re = joblib.load(os.path.join(models_dir, 'reg_energy_' + str(telescope_model) + '.sav'))
        logging.info('Loading model: ' + models_dir + '/reg_energy_' + str(telescope_model) + '.sav')
        logging.info('Energy reconstruction..')
        dl2['log_reco_energy'] = energy_re.predict(dl1[energy_regression_features])
        dl1['log_reco_energy'] = dl2['log_reco_energy']
        dl2['reco_energy'] = 10**dl2['log_reco_energy']
    except FileNotFoundError:
        logging.warning('Energy reconstruction model was not found. Energy of all reconstructed events will be set to 0.')
        dl2['log_reco_energy'] = 0
        dl1['log_reco_energy'] = 0
        dl2['reco_energy'] = 0

    if (config['disp_method'] == 'disp_vector') and not stereo:

        logging.info('Mono DISP reconstruction using disp_vector..')

        try:
            logging.info('Loading model: ' + models_dir + '/reg_disp_' + str(telescope_model) + '.sav')
            disp_re = joblib.load(os.path.join(models_dir, 'reg_disp_' + str(telescope_model) + '.sav'))
            disp_vec = disp_re.predict(dl1[disp_regression_features])
        except FileNotFoundError:
            logging.warning('DISP model was not found. DISP vector of all reconstructed events will be set to [0, 0].')
            disp_vec = np.zeros((len(dl2), 2))


    elif (config['disp_method'] == 'disp_norm_sign') and not stereo:

        logging.info('Mono DISP reconstruction using disp_norm_sign..')
        try:
            logging.info('Loading model: ' + models_dir + '/reg_disp_norm_' + str(telescope_model) + '.sav')
            disp_re = joblib.load(os.path.join(models_dir, 'reg_disp_norm_' + str(telescope_model) + '.sav'))
            logging.info('Loading model: ' + models_dir + '/class_disp_sign_' + str(telescope_model) + '.sav')
            disp_clas = joblib.load(os.path.join(models_dir, 'class_disp_sign_' + str(telescope_model) + '.sav'))

            disp_norm = disp_re.predict(dl1[disp_regression_features])
            disp_sign = disp_clas.predict(dl1[disp_classification_features])

            dl2['reco_disp_norm'] = disp_norm
            dl1['reco_disp_norm'] = disp_norm
            dl2['reco_disp_sign'] = disp_sign

            # NOTE: Be carefull with the units! - must be revisited for proper unit conversion!
            disp_angle = dl2['camera_frame_hillas_psi'] * np.pi / 180

            disp_vec = disp_vector(disp_norm, disp_angle, disp_sign)

        except FileNotFoundError:
            logging.warning('DISP model was not found. DISP vector of all reconstructed events will be set to [0, 0].')
            disp_vec = np.zeros((len(dl2), 2))

    elif stereo:

        logging.info('Mono DISP reconstruction of disp norm..')
        disp_re = joblib.load(os.path.join(models_dir, 'reg_disp_norm_' + str(telescope_model) + '.sav'))
        disp_norm = disp_re.predict(dl1[disp_regression_features])
        dl2['reco_disp_norm'] = disp_norm
        dl1['reco_disp_norm'] = disp_norm

        # For MARS-like stereo reconstruction we need to consider both possible source positions
        dl2['reco_disp_sign_p'] = 1
        dl2['reco_disp_sign_m'] = -1

    else:
        raise Exception("ERROR: disp_method in the config file not recognized, use one of these: disp_vector, disp_norm_sign")

    # g/h needs to be reconstructed here, after reco_log_energy and reco_disp_norm columns are in dl1
    try:
        logging.info('Loading model: ' + models_dir + '/class_gh_' + str(telescope_model) + '.sav')
        gh_clas = joblib.load(os.path.join(models_dir, 'class_gh_' + str(telescope_model) + '.sav'))
        logging.info('Gamma/hadron classification..')
        dl2['gammaness'] = gh_clas.predict_proba(dl1[gh_clas_features])[:, 0]
    except FileNotFoundError:
        logging.warning('Gamma/hadron classification model was not found. Gammaness of all reconstructed events will be set to 0.')
        dl2['gammaness'] = 0

    if not stereo:

        dl2['reco_disp_dx'] = disp_vec[:, 0]
        dl2['reco_disp_dy'] = disp_vec[:, 1]
        dl2['reco_src_x'], dl2['reco_src_y'] = disp_to_pos(dl2.reco_disp_dx,
                                                            dl2.reco_disp_dy,
                                                            dl2.camera_frame_hillas_x,
                                                            dl2.camera_frame_hillas_y,
                                                            )

        # NOTE: Be carefull with the units! - must be revisited for proper unit conversion!
        if ~(dl2.true_alt_tel.values > np.pi).any():
            raise Exception("ERROR: true_alt_tel, true_az_tel are probably not in degrees! Check your DL1 file.")

        if mc:
            times=None
        else:
            times = Time(dl2.local_time, format='unix', scale='utc')

        src_pos_reco = reco_source_position_sky(dl2.camera_frame_hillas_x.values * u.m,
                                                    dl2.camera_frame_hillas_y.values * u.m,
                                                    dl2.reco_disp_dx.values * u.m,
                                                    dl2.reco_disp_dy.values * u.m,
                                                    dl2.equivalent_focal_length.values * u.m,
                                                    dl2.true_alt_tel.values * u.deg,
                                                    dl2.true_az_tel.values * u.deg,
                                                    location=telescope_location(subarray, telescope),
                                                    times=times)

        dl2['reco_alt'] = src_pos_reco.alt.deg
        azimuth_transformed = (src_pos_reco.az.to_value(u.deg) +  180) % 360 - 180
        if abs(dl2['true_az_tel'][0]) <= 90:
            dl2['reco_az'] = azimuth_transformed
            logging.info('Using transformed azimuth')
        else:
            dl2['reco_az'] = src_pos_reco.az.deg
            logging.info('Using non-transformed azimuth')

        radec = src_pos_reco.transform_to('icrs')
        dl2['reco_ra'] = radec.ra.deg
        dl2['reco_dec'] = radec.dec.deg

        logging.info('Reconstruction done.')

    return dl2


def stereo_reconstruction(
        params, config=None,
        ismc=False, telescopes=None, subarray=None):
    """
    Averages energy and gammaness stereoscopicly reconstructed for
    each telescope to get final stereo quantities for each event.
    For direction reconstruction, alt/az for both signs from each
    telescope are calculated. Then the combination leading to
    reconstructed coordinates with the smallest distance is
    selected and final direction is calculated as weighted average
    of the two (MARS-like stereo reconstruction).

    Parameters
    ----------
    params: pandas.DataFrame
    config: dict
    ismc: bool
    telescopes: list of strings
    subarray: ctapipe.instrument.SubarrayDescription
        subarray of the input file (``SubarrayDescription.from_hdf(input_file)``), for the telescope positions

    Returns
    -------
    pandas.DataFrame

    """
    if subarray is None:
        raise ValueError("subarray is required (SubarrayDescription.from_hdf(input_file)): the telescope positions")

    logging.info('Stereo reconstruction..')
    weights = config["analysis"]["stereo_reco_weights"]
    logging.info('Stereo reconstruction event weights: %s', weights)

    energy_average, energy_var, energy_wvar = get_average_param(params, param='log_reco_energy', weights=weights)
    gammaness_average, gammaness_var, gammaness_wvar = get_average_param(params, param='gammaness', weights=weights)

    dl2 = get_stereo_dl2(params, ismc=ismc)

    dl2['log_reco_energy'] = energy_average
    dl2['reco_energy'] = 10**energy_average
    dl2['var_reco_energy'] = (np.log(10)*10**energy_average)**2 * energy_var
    dl2['wvar_reco_energy'] = (np.log(10)*10**energy_average)**2 * energy_wvar
    dl2['gammaness'] = gammaness_average
    dl2['var_gammaness'] = gammaness_var
    dl2['wvar_gammaness'] = gammaness_wvar

    # We need to add these for event types classificantion
    if ismc:
        mask_tel1 = params['tel_id'] == 1
        mask_tel2 = params['tel_id'] == 2
    else:
        mask_tel1 = params['tel_id'] == 21
        mask_tel2 = params['tel_id'] == 22
    hillas_tel1 = params[mask_tel1][["obs_id", "event_id", "HillasReconstructor_tel_impact_distance"]]
    hillas_tel2 = params[mask_tel2][["obs_id", "event_id", "HillasReconstructor_tel_impact_distance"]]
    hillas_tel1 = hillas_tel1.rename(
        columns={"HillasReconstructor_tel_impact_distance":
        "HillasReconstructor_tel_impact_distance_tel1"}
        )
    hillas_tel2 = hillas_tel2.rename(
        columns={"HillasReconstructor_tel_impact_distance":
        "HillasReconstructor_tel_impact_distance_tel2"}
        )
    dl2 = dl2.merge(hillas_tel1, on=['obs_id', 'event_id'], how='left')
    dl2 = dl2.merge(hillas_tel2, on=['obs_id', 'event_id'], how='left')

    # Arrival direction reconstruction

    # NOTE: Be carefull with the units! - must be revisited for proper unit conversion!
    disp_angle = params['camera_frame_hillas_psi'] * np.pi / 180

    # Here we reconstruct both possible alt,az of the event for each telescope
    for sign in ['reco_disp_sign_p', 'reco_disp_sign_m']:

        disp_sign = params[sign]

        # initialization of two more columns
        params['reco_alt'+sign.replace('reco_disp','')] = 0
        params['reco_az'+sign.replace('reco_disp','')] = 0

        disp_vec = disp_vector(params['reco_disp_norm'], disp_angle, disp_sign)
        reco_disp_dx = disp_vec[:, 0]
        reco_disp_dy = disp_vec[:, 1]

        #for telescope in ['tel_001', 'tel_002']:
        # This should hopefuly work on data as well (where we have tel_021, tel_022)
        for telescope in telescopes:

            tel = get_tel_id(telescope)
            mask = params['tel_id'] == tel

            src_pos_reco = reco_source_position_sky(
                params[mask].camera_frame_hillas_x.values * u.m,
                params[mask].camera_frame_hillas_y.values * u.m,
                reco_disp_dx[mask] * u.m,
                reco_disp_dy[mask] * u.m,
                params[mask].equivalent_focal_length.values * u.m,
                params[mask].true_alt_tel.values * u.deg,
                params[mask].true_az_tel.values * u.deg,
                location=subarray.tel_earth_locations[tel],
                )

            azimuth_transformed = (src_pos_reco.az.to_value(u.deg) + 180) % 360 - 180
            params.loc[mask, 'reco_alt'+sign.replace('reco_disp','')] = src_pos_reco.alt.deg
            if abs(params[mask].true_az_tel.values[0]) <= 90:
                params.loc[mask, 'reco_az'+sign.replace('reco_disp','')] = azimuth_transformed
                logging.info('Using transformed azimuth')
            else:
                params.loc[mask, 'reco_az'+sign.replace('reco_disp','')] = src_pos_reco.az.deg
                logging.info('Using non-transformed azimuth')

    # Now we have in params two azimuths and altitudes for each telescope, generated
    # using either + or - sign togther with the reconstructd disp_norm. We now need
    # to calculate all four possible angular distances (+,+), (-,-), (+,-), (-,+)
    # select the closest one, and calculate weighted average direction
    # NOTE: We can safely average directions in alt/az, no need to go for ra/dec, because the difference
    # between the two frames is minor (of the order of 10^-3 deg).

    averaged_direction = get_averaged_direction(params, weights=weights, telescopes=telescopes)

    dl2 = pd.merge(dl2, averaged_direction, on=["obs_id", "event_id"])

    # Here we can use horizon frame of either of the telescopes. The simulations have no times:
    # an arbitrary time is used
    location = telescope_location(subarray, telescopes[0])
    if ismc:
        horizon_frame = AltAz(location=location, obstime=SIMULATION_OBSTIME)
    else:
        horizon_frame = AltAz(location=location, obstime=Time(dl2['local_time'], format='unix'))
    dl2 = add_reco_ra_dec(dl2, horizon_frame=horizon_frame)

    tel = []
    for tele in telescopes:
        tel.append(int(tele.split('_')[-1]))

    dl2.drop([
        'reco_alt_sign_p_tel'+str(tel[0]),
        'reco_alt_sign_m_tel'+str(tel[0]),
        'reco_az_sign_p_tel'+str(tel[0]),
        'reco_az_sign_m_tel'+str(tel[0]),
        'reco_alt_sign_p_tel'+str(tel[1]),
        'reco_alt_sign_m_tel'+str(tel[1]),
        'reco_az_sign_p_tel'+str(tel[1]),
        'reco_az_sign_m_tel'+str(tel[1])
        ],
        inplace=True,
        axis=1
        )

    return dl2


def get_averaged_direction(params, weights=None, telescopes=None):
    """
    Calculate weighted average of the closest from all four possible
    alt/az combinations from both telescopes.

    Parameters
    ----------
    params: pandas.DataFrame
    weights: str
        Usualy one of the hillas features present in the input
        DL1 pandas.DataFrame
    telescopes: list of strings

    Returns
    -------
    pandas.DataFrame

    """

    tel = []
    for tele in telescopes:
        tel.append(int(tele.split('_')[-1]))

    params_reshaped = reshape_event_table_for_stereo(params, tels=tel)

    logging.info('Finding the closest combination of direction..')
    distance = pd.DataFrame()
    for combo in ['pp', 'mm', 'mp', 'pm']:
        distance[combo] = angular_distance(params_reshaped, combination=combo, tels=tel)

    params_reshaped['min_distance'] = distance.min(axis = 1)
    params_reshaped['best_sign_combination'] = distance.idxmin(axis = 1)

    params_reshaped['reco_alt'] = 0
    params_reshaped['reco_az'] = 0

    telescopes_averaged = params_reshaped.copy(deep=True)

    logging.info('Averaging reconstructed event directions..')
    for combo in ['pp', 'mm', 'mp', 'pm']:
        mask = telescopes_averaged['best_sign_combination'] == combo
        for coord in ['alt', 'az']:
            if weights is not None:
                telescopes_averaged.loc[mask, 'reco_'+coord] = (
                    params_reshaped.loc[mask]['reco_'+coord+'_sign_'+combo[0]+'_tel'+str(tel[0])] *
                    params_reshaped.loc[mask][weights + '_tel'+str(tel[0])] +
                    params_reshaped.loc[mask]['reco_'+coord+'_sign_'+combo[1]+'_tel'+str(tel[1])] *
                    params_reshaped.loc[mask][weights + '_tel'+str(tel[1])]
                    ) / (
                        params_reshaped.loc[mask][weights + '_tel'+str(tel[0])] +
                        params_reshaped.loc[mask][weights + '_tel'+str(tel[1])]
                    )
            else:
                telescopes_averaged.loc[mask, 'reco_'+coord] = (
                    params_reshaped.loc[mask]['reco_'+coord+'_sign_'+combo[0]+'_tel'+str(tel[0])] +
                    params_reshaped.loc[mask]['reco_'+coord+'_sign_'+combo[1]+'_tel'+str(tel[1])]
                    ) / 2

    return telescopes_averaged


def reshape_event_table_for_stereo(params, tels=None):

    data_tel1 = get_data_tel(params, tel=tels[0])
    data_tel2 = get_data_tel(params, tel=tels[1])
    params_reshaped = pd.merge(data_tel1, data_tel2, on=["obs_id", "event_id"])

    return params_reshaped


def get_data_tel(params, tel=1):

    mask_tel = params['tel_id'] == tel

    # We need those telescope dependent parameters (e.g. length or width)
    # for later use for event selection in performance evaluation
    data_tel = params[mask_tel][[
        "obs_id",
        "event_id",
        'reco_alt_sign_p',
        'reco_alt_sign_m',
        'reco_az_sign_p',
        'reco_az_sign_m',
        "camera_frame_hillas_intensity",
        "camera_frame_hillas_width",
        "camera_frame_hillas_length",
        "leakage_intensity_width_2",
        "camera_frame_hillas_wl",
        "camera_frame_hillas_skewness",
        "camera_frame_hillas_kurtosis",
        "camera_frame_timing_slope"
        ]]

    data_tel = data_tel.rename(
                    columns={
                            "reco_alt_sign_p": "reco_alt_sign_p_tel"+str(tel),
                            "reco_alt_sign_m": "reco_alt_sign_m_tel"+str(tel),
                            "reco_az_sign_p": "reco_az_sign_p_tel"+str(tel),
                            "reco_az_sign_m": "reco_az_sign_m_tel"+str(tel),
                            "camera_frame_hillas_intensity": "camera_frame_hillas_intensity_tel"+str(tel),
                            "camera_frame_hillas_width": "camera_frame_hillas_width_tel"+str(tel),
                            "camera_frame_hillas_length": "camera_frame_hillas_length_tel"+str(tel),
                            "leakage_intensity_width_2": "leakage_intensity_width_2_tel"+str(tel),
                            "camera_frame_hillas_wl": "camera_frame_hillas_wl_tel"+str(tel),
                            "camera_frame_hillas_skewness": "camera_frame_hillas_skewness_tel"+str(tel),
                            "camera_frame_hillas_kurtosis": "camera_frame_hillas_kurtosis_tel"+str(tel),
                            "camera_frame_timing_slope": "camera_frame_timing_slope_tel"+str(tel)
                            }
                )

    return data_tel


def angular_distance(params, combination='pp', tels=None):
    """
    Calculates angular distance for all four combination of
    reconstructed source positions for images of the two telescopes.

    Parameters
    ----------
    params: pandas.DataFrame
    combination: str
    tels: list of ints [1,2]

    Returns
    -------
    pandas.DataFrame
        distance in radians
    """

    distance = angular_separation_altaz(
                    params['reco_alt_sign_'+combination[0]+'_tel'+str(tels[0])].values * np.pi/180. * u.rad,
                    params['reco_az_sign_'+combination[0]+'_tel'+str(tels[0])].values * np.pi/180. * u.rad,
                    params['reco_alt_sign_'+combination[1]+'_tel'+str(tels[1])].values * np.pi/180. * u.rad,
                    params['reco_az_sign_'+combination[1]+'_tel'+str(tels[1])].values * np.pi/180. * u.rad,
                )
    return distance


def get_average_param(params, param=None, weights=None):
    """
    Get stereo quantities for both telescopes averaged

    Parameters
    ----------
    params: pandas.DataFrame
    param: str
        Quantity to be averaged, e.g. gammaness
    weights: str

    Returns
    -------
    average:
        pandas.DataFrame
    variance:
        Note that variance is calculated always from unweighted quantities
        pandas.DataFrame
    """

    logging.info('Averaging reconstructed %s', param)

    if weights is not None:
        weights = params[weights]

        df = pd.DataFrame(
            data={
                "obs_id": params["obs_id"],
                "event_id": params["event_id"],
                "weight": weights,
                "param": params[param],
                "weighted_param": params[param] * weights,
            }
        )

        df_sum = df.groupby(["obs_id", "event_id"]).sum()
        average = df_sum["weighted_param"] / df_sum["weight"]
        variance = df.groupby(["obs_id", "event_id"]).var()["param"]

        # weighted variance
        # Merge the mean back to the original df for subtraction
        df = df.join(average.rename("mean_weighted"), on=["obs_id", "event_id"])
        # Calculate (x - mean)^2 * weight
        df["squared_diff"] = df["weight"] * (df["param"] - df["mean_weighted"])**2
        # Sum over groups and divide by sum of weights to get variance
        sum_w_squared_diff = df.groupby(["obs_id", "event_id"])["squared_diff"].sum()
        weighted_variance = sum_w_squared_diff / df_sum["weight"]

    else:
        df = pd.DataFrame(
            data={
                "obs_id": params["obs_id"],
                "event_id": params["event_id"],
                "param": params[param],
            }
        )

        df_sum = df.groupby(["obs_id", "event_id"]).sum()
        average = df_sum["param"] / df.groupby(["obs_id", "event_id"]).size()
        variance = df.groupby(["obs_id", "event_id"]).var()["param"]
        variance = weighted_variance

    return average, variance, weighted_variance


def get_stereo_dl2(params, ismc=False):
    """
    Select columns and create DL2 stereo table

    Parameters
    ----------
    params: pandas.DataFrame
    ismc: bool

    Returns
    -------
    pandas.DataFrame

    """

    features = [
        "obs_id",
        "event_id",
        "true_az_tel",
        "true_alt_tel",
        "HillasReconstructor_core_x",
        "HillasReconstructor_core_y",
        "HillasReconstructor_h_max"
        ]

    # For MC we store some extra true parameters in DL2 for performance evaluation
    if ismc:
        features = features + [
            "true_az",
            "true_alt",
            "true_energy",
            "log_true_energy",
            "true_core_x",
            "true_core_y",
            "true_h_first_int",
            "true_x_max",
            "true_shower_primary_id",
            "true_camera_x",
            "true_camera_y",
        ]
    else:
        # For data we need to add also some more parameters
        features = features + [
            "local_time",
            "date",
            "equivalent_focal_length"
            ]

    if 'min_true_energy_cut' in params.keys():
        features.append('min_true_energy_cut')

    # Here we ask for mean(), but these quantities should be the same for both telescopes
    # NOTE: ..Maybe except true_camera_x, y? This should be checked!
    dl2 = params[features].groupby(["obs_id", "event_id"]).mean()
    return dl2
