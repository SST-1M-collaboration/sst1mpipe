"""
Daily analysis of the data of a night on Calculus, per target: R0 -> DL1 (``sst1mpipe-monitoring``,
``sst1mpipe-process``) for each telescope, DL1 stereo (``sst1mpipe_data_dl1_dl1_stereo``), DL1 -> DL2
(``sst1mpipe-apply-models``, models of the closest node of ``mono_model_dir`` / ``stereo_model_dir``).

TODO: DL2 -> DL3 and the DL1 distributions from the ctapipe DL2 (``sst1mpipe_data_dl2_dl3`` and
``sst1mpipe_extract_dl1_distributions`` read the DL2 of the former sst1mpipe random forests).

Usage::

    nohup python daily_analysis.py /path/to/daily_ana_config.json 2>daily_ana.err 1>daily_ana.out &
    python daily_analysis.py  # last night, default configuration
"""
import datetime
import glob
import json
import logging
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
from functools import partial
from pathlib import Path

from astropy.io import fits

from sst1mpipe.resources import MONITORING_CONFIG_FILE

DEFAULT_CONFIG_FILE = '/data/work/analysis/Daily_analysis/default_daily_ana_config.json'
IGNORED_TARGETS = {'DARK', 'Transition', 'UNKNOWN', 'BIAS', 'WRtest', 'TRANSITION', 'transition', 'dark'}


def run(cmd):
    logging.info(' '.join(map(str, cmd)))
    subprocess.run(list(map(str, cmd)), check=False)


def r0_dl1(input_file, config_file, monitoring_config_file, out_dir):
    """R0 -> DL1 of a raw file: monitoring, then calibration and DL1"""
    run_name = Path(input_file).name.removesuffix('.fits.fz')
    monitoring = Path(out_dir) / f'{run_name}.monitoring.h5'
    run(['sst1mpipe-monitoring', '-i', input_file, '-o', monitoring, '--config', monitoring_config_file,
         '--overwrite'])
    run(['sst1mpipe-process', '-i', input_file, '-o', Path(out_dir) / f'{run_name}.dl1.h5', '--config', config_file,
         '--monitoring-source', 'SST1MMonitoringSource', '--SST1MMonitoringSource.input_files', monitoring,
         '--overwrite', '--no-progress'])


def dl1_dl1_stereo(input_file, config_file, tel2_dir, out_dir):
    """coincident events of a DL1 file of tel 1 with the DL1 files of tel 2"""
    run(['sst1mpipe_data_dl1_dl1_stereo', '--input-file-tel1', input_file, '--config', config_file,
         '--input-dir-tel2', tel2_dir, '--output-dir', out_dir])


def dl1_dl2(input_file, models_dir, out_dir):
    """DL1 -> DL2 (ctapipe DL2) with the models of the node closest to the pointing"""
    output = Path(out_dir) / Path(input_file).name.replace('dl1', 'dl2')
    run(['sst1mpipe-apply-models', '-i', input_file, '-m', models_dir, '-o', output, '--overwrite',
         '--no-progress'])


def pool_map(n_proc, function, files, **kwargs):
    with mp.Pool(n_proc) as pool:
        pool.map(partial(function, **kwargs), files)


def make_runlist_allfiles(itel,year,month,day,rootdir='/net/'):
        """
        make a list of raw data file path

        Parameters
        ----------
        itel : int
            telescope number (1 or 2)
        year  : int

        month : int

        day   : int

        Returns
        -------
        list of raw data file paths
        """
        basedir = rootdir+f'/cs{itel}/data/raw/'
        filerad = f'SST1M{itel}'
        filedir = os.path.join(basedir,
                                       str(year),
                                       f"{month:02d}",
                                       f"{day:02d}",
                                       filerad
                                       )
        return glob.glob(filedir+"/"+filerad+"*.fits.fz")

def refine_file_list(file_list):
    """
       sort list of raw data fits based on there TARGET field in header

        Parameters
        ----------
        file_list : list
            list of raw data file paths

        Returns
        -------
        dict of raw data file paths
    """

    list_dict = dict()
    list_dict['DARK']       = []
    list_dict['Transition'] = []
    list_dict['UNKNOWN']    = []
    for fitsfile in file_list:
        f = fits.open(fitsfile)
        try:
            target=f[2].header['TARGET'].split(',')[0].split('_')[0].replace(" ", "")
        except Exception:
            target='UNKNOWN'
        if target not in list_dict.keys():
            list_dict[target]    = []
        list_dict[target].append(fitsfile)
        f.close()
    return list_dict

def run_daily_ana(daily_config):
    """R0 -> DL2 of the data of a night (``daily_config``: json configuration)"""
    n_proc = daily_config["n_proc"]
    ana_dir = daily_config["ana_dir"]
    config_file = daily_config["config_file"]
    monitoring_config_file = daily_config.get("monitoring_config_file", MONITORING_CONFIG_FILE)
    year, month, day = daily_config.get("year"), daily_config.get("month"), daily_config.get("day")
    model_dirs = {"cs1": daily_config["mono_model_dir"], "cs2": daily_config["mono_model_dir"],
                  "stereo": daily_config["stereo_model_dir"]}
    stereo = daily_config["stereo"]
    mono = {"cs1": daily_config["cs1"], "cs2": daily_config["cs2"]}
    run_dl1, run_dl2 = daily_config["run_dl1"], daily_config["run_dl2"]

    if (year is None) or (month is None) or (day is None):
        yesterday = datetime.datetime.now() - datetime.timedelta(days=1)
        year, month, day = yesterday.year, yesterday.month, yesterday.day
    date = f'{year:04d}{month:02d}{day:02d}'

    Path(ana_dir, 'logs').mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(Path(ana_dir, 'logs', f'daily_ana_{date}.log'), 'w+'),
                  logging.StreamHandler(stream=sys.stdout)],
    )
    logging.info("Daily analysis (%s) START", date)
    if any(daily_config.get(key) for key in ("run_dl3", "run_dl3_eff", "run_data_qual", "run_data_qual_eff")):
        logging.warning("DL3 and DL1 distributions not available yet from the ctapipe DL2: skipped")

    files = {tel: refine_file_list(make_runlist_allfiles(itel, year, month, day, rootdir=daily_config["root_dir"]))
             for itel, tel in ((1, "cs1"), (2, "cs2"))}
    if not any(len(v) for d in files.values() for v in d.values()):
        logging.warning("NO DATA FOUND : daily analysis ended")
        return
    if not all(any(len(v) for v in files[tel].values()) for tel in ("cs1", "cs2")):
        stereo = False

    datedir = Path(ana_dir, date)
    datedir.mkdir(exist_ok=True)
    with open(datedir / 'daily_config.json', 'w') as f:
        json.dump(daily_config, f)
    shutil.copy(config_file, datedir)

    for target in sorted((set(files["cs1"]) | set(files["cs2"])) - IGNORED_TARGETS):
        logging.info("%s: CS1 %d files, CS2 %d files", target, len(files["cs1"].get(target, [])),
                     len(files["cs2"].get(target, [])))
        # DL1 of each telescope for its mono analysis or the stereo, DL2 of the analyses asked
        dl2 = {**mono, "stereo": stereo}
        dirs = {}
        for tel in [tel for tel in ("cs1", "cs2") if mono[tel] or stereo] + (["stereo"] if stereo else []):
            dirs[tel] = {level: datedir / target / tel / level for level in ("dl1", "dl2")}
            for path in dirs[tel].values():
                path.mkdir(parents=True, exist_ok=True)

        for tel in ("cs1", "cs2"):
            if tel in dirs and run_dl1 and files[tel].get(target):
                pool_map(n_proc, r0_dl1, files[tel][target], config_file=config_file,
                         monitoring_config_file=monitoring_config_file, out_dir=dirs[tel]["dl1"])
        if stereo and run_dl1:
            pool_map(n_proc, dl1_dl1_stereo, glob.glob(str(dirs["cs1"]["dl1"] / '*.dl1.h5')),
                     config_file=config_file, tel2_dir=dirs["cs2"]["dl1"], out_dir=dirs["stereo"]["dl1"])

        if run_dl2:
            for tel in [tel for tel in dirs if dl2[tel]]:
                pattern = '*.h5' if tel == "stereo" else '*.dl1.h5'  # not the monitoring files
                pool_map(n_proc, dl1_dl2, glob.glob(str(dirs[tel]["dl1"] / pattern)),
                         models_dir=model_dirs[tel], out_dir=dirs[tel]["dl2"])

    logging.info("Daily analysis (%s) ended", date)


if __name__ == '__main__':
    with open(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CONFIG_FILE) as json_data_file:
        run_daily_ana(json.load(json_data_file))
