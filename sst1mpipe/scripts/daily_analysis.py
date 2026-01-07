import os
import multiprocessing as mp
import numpy as np
import glob
import copy
import sys
import json
import shutil
import logging
import datetime 

from pathlib import Path
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

'''
this run the analysis pipeline up to stereo dl2 on Calculus.
USAGE :

nohup python daily_analysis.py /path/to/some/ana_config.json 2>daily_ana.err 1>daily_ana.out &

or 

for auto proccessing of the last night data :
python daily_analysis.py 2>daily_ana.err 1>daily_ana.out

'''


DEFAULT_CONFIG_FILE = '/data/work/analysis/Daily_analysis/default_daily_ana_config.json'

class iargs:
    """
    dummy object made to handle arguments to 
    be passed to the multiprocessing pool.
    """
    def __init__(self):
        pass
        

def args_maker(arg,file_list):
    '''
    Create a list of dummy arg object to be passed to the multiprocessing pool.
        Parameters
        ----------
        arg: class iargs obj.
            dummy arg object that contain arguments to be
            padded to sst1mpipe scripts

        file_list: list
            list of file paths 
    '''

    args_list = []
    
    for input_file in file_list:
        t_arg = copy.copy(arg)
        t_arg.input_file=input_file
        args_list.append(t_arg)
    return args_list

    
    

def r0_dl1_1file(arg):
        """
        python binding of sst1mpipe_r0_dl1
        """

        cmd = 'sst1mpipe_r0_dl1 --input-file {} \
               --config {} \
               --output-dir {} --precise-timestamps'.format(arg.input_file,
                                                            arg.config_file,
                                                            arg.out_dir)
        logging.info(cmd)
        os.system(cmd)

def dl1_dl1_1file(arg):
        """
        python binding of sst1mpipe_data_dl1_dl1_stereo
        """
        cmd = 'sst1mpipe_data_dl1_dl1_stereo \
               --input-file {} \
               --config {} \
               --input-dir-tel2 {} \
               --output-dir {}'.format(arg.input_file,
                                       arg.config_file,
                                       arg.tel2_dir,
                                       arg.out_dir)
        logging.info(cmd)
        os.system(cmd)


def dl1_dl2_1file(arg):
        """
        python binding of sst1mpipe_dl1_dl2
        """
        cmd = 'sst1mpipe_dl1_dl2 \
               --input-file {} \
               --config {} \
               --models-dir {} \
               --output-dir {} --stereo'.format(arg.input_file,
                                         arg.config_file,
                                         arg.models_dir,
                                         arg.out_dir)
        logging.info(cmd)
        os.system(cmd)


def dl1_dl2_1file_1tel(arg):
        """
        python binding of sst1mpipe_dl1_dl2
        """

        cmd = 'sst1mpipe_dl1_dl2 \
               --input-file {} \
               --config {} \
               --models-dir {} \
               --output-dir {}'.format(arg.input_file,
                                       arg.config_file,
                                       arg.models_dir,
                                       arg.out_dir)
        logging.info(cmd)
        os.system(cmd)

def dl2_dl3_1dir(arg):
        """
        python binding of sst1mpipe_data_dl2_dl3
        """

        cmd = 'sst1mpipe_data_dl2_dl3 \
               --input-dir {} \
               --irf-dir {} \
               --config {} \
               --target-name {} \
               --target-ra {} \
               --target-dec {} \
               --output-dir {}'.format(arg.input_dir,
                                       arg.irf_dir,
                                       arg.config_file,
                                       arg.target_name,
                                       arg.target_ra,
                                       arg.target_dec,
                                       arg.out_dir)
        if arg.gammaness_cut_dir is not None:
            cmd = cmd + " --gammaness-cut-dir {}".format(arg.gammaness_cut_dir)
        logging.info(cmd)
        os.system(cmd)


def extract_dl1_distributions(arg):
        """
        python binding of sst1mpipe_extract_dl1_distributions
        """

        cmd = 'sst1mpipe_extract_dl1_distributions \
               --dl1-dir {} \
               --date {} \
               --output-dir {} \
               --histogram-bins 100 \
               --dl3-index-dir {}'.format(arg.dl1_dir,
                                          arg.date_str,
                                          arg.out_dir,
                                          arg.dl3_dir)
        logging.info(cmd)
        os.system(cmd)


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
        basedir = rootdir+'/cs{}/data/raw/'.format(itel)
        filerad = 'SST1M{}'.format(itel)
        filedir = os.path.join(basedir,
                                       str(year),
                                       "{:02d}".format(month),
                                       "{:02d}".format(day),
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
    """
    run daily r0 to dl3 analysis 

    Parameters
    ----------
    daily_config : json config object
            
    Returns
    -------
    None
    """

    n_proc             = daily_config["n_proc"]
    ana_dir            = daily_config["ana_dir"]
    config_file        = daily_config["config_file"]
    year               = daily_config["year"]
    month              = daily_config["month"]
    day                = daily_config["day"]
    mono_model_dir     = daily_config["mono_model_dir"]
    stereo_model_dir   = daily_config["stereo_model_dir"]
    cs1                = daily_config["cs1"]
    cs2                = daily_config["cs2"]
    stereo             = daily_config["stereo"]
    rootdir            = daily_config["root_dir"]
    run_dl1            = daily_config["run_dl1"]
    run_dl2            = daily_config["run_dl2"]
    run_dl3            = daily_config["run_dl3"]
    run_data_qual      = daily_config["run_data_qual"]

    run_dl3_eff        = daily_config["run_dl3_eff"]
    run_data_qual_eff  = daily_config["run_data_qual_eff"]
    eff_cut_dir_mono   = daily_config["eff_cut_dir_mono"]
    eff_cut_dir_stereo = daily_config["eff_cut_dir_stereo"]
    try:
        output_logfile = os.path.join(ana_dir, "logs", 'daily_ana_{:04d}{:02d}{:02d}.log'.format(year,month,day))
    except Exception:
        output_logfile = os.path.join(ana_dir, "logs", 'default_daily_ana.log')
    
    Path(ana_dir + '/logs/').mkdir(exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers= [
            logging.FileHandler(output_logfile, 'w+'),
            logging.StreamHandler(stream=sys.stdout)
            ]
    )
    
    if (year is None) or (month is None) or (day is None):
        logging.info("Running daily analysis for last night !")
        year  = (datetime.datetime.now()-datetime.timedelta(days=1)).year
        month = (datetime.datetime.now()-datetime.timedelta(days=1)).month
        day   = (datetime.datetime.now()-datetime.timedelta(days=1)).day
    logging.info("---------------------------------------------")
    logging.info("Daily analysis ({:04d}/{:02d}/{:02d}) START ".format(year,month,day))
    
    ## get raw file list 
    raw_file_list_t1 = make_runlist_allfiles(itel = 1,
                                             year=year,
                                             month=month,
                                             day=day,
                                             rootdir=rootdir)

    raw_file_list_t2 = make_runlist_allfiles(itel = 2,
                                             year=year,
                                             month=month,
                                             day=day,
                                             rootdir=rootdir)
 
    
    if (len(raw_file_list_t1)==0) and (len(raw_file_list_t2)==0):
        logging.warning("NO DATA FOUND : daily analysis ended")
        return

    if (len(raw_file_list_t1)==0) or (len(raw_file_list_t2)==0):
        stereo = False



    ## refine file list by target in header
    dict_list_t1 = refine_file_list(raw_file_list_t1)
    dict_list_t2 = refine_file_list(raw_file_list_t2)
    
    
    target_list = np.unique(list(dict_list_t1.keys())+list(dict_list_t2.keys()))
    logging.info("Target found ::")
    for target in target_list:
        logging.info('{} :'.format(target))
        try :
            logging.info('CS1 : {}'.format(len(dict_list_t1[target]) ))
        except Exception:
           logging.warning('no CS1 data')
        try :
            logging.info('CS2 : {}'.format(len(dict_list_t2[target]) ))
        except Exception:
           logging.warning('no CS2 data')

 

    ## make direcotory
    datedir = os.path.join(ana_dir,
                           '{:04d}{:02d}{:02d}'.format(year,month,day) )
    Path(datedir).mkdir(exist_ok=True)

    ## copy config files
    with open(datedir+'/daily_config.json', 'w') as f:
        json.dump(daily_config, f)

    shutil.copy(config_file, datedir)

     
    for target in target_list:
        if target in ['DARK','Transition','UNKNOWN','BIAS','WRtest',"TRANSITION","transition","dark"]:
            continue
        target_dir = os.path.join(datedir,
                                  '{}'.format(target) )
        Path(target_dir).mkdir(exist_ok=True)

        ## yes, this have a low probability to work!
        ## TODO get target pos.. How ? 
        ## is ok, we don't realy use this for now
        try:
            target_pos = SkyCoord.from_name(target)
        except Exception:
            logging.info("Cant guess the target position!")
            target_pos = SkyCoord(0*u.deg,0*u.deg)

        #CS1
        if cs1 or stereo:
            cs1_dir = os.path.join(target_dir,'cs1')
            Path(cs1_dir).mkdir(exist_ok=True)

            cs1_dl1_dir = os.path.join(cs1_dir,'dl1')
            Path(cs1_dl1_dir).mkdir(exist_ok=True)

            cs1_dl2_dir = os.path.join(cs1_dir,'dl2')
            Path(cs1_dl2_dir).mkdir(exist_ok=True)

            cs1_dl3_dir = os.path.join(cs1_dir,'dl3')
            Path(cs1_dl3_dir).mkdir(exist_ok=True)

            cs1_dqual_dir = os.path.join(cs1_dir,'distributions')
            Path(cs1_dqual_dir).mkdir(exist_ok=True)

            cs1_dl3_dir_eff = os.path.join(cs1_dir,'dl3_eff')
            Path(cs1_dl3_dir_eff).mkdir(exist_ok=True)

            cs1_dqual_dir_eff = os.path.join(cs1_dir,'distributions_eff')
            Path(cs1_dqual_dir_eff).mkdir(exist_ok=True)

        ## CS2
        if cs2 or stereo:
            cs2_dir = os.path.join(target_dir,'cs2')
            Path(cs2_dir).mkdir(exist_ok=True)

            cs2_dl1_dir = os.path.join(cs2_dir,'dl1')
            Path(cs2_dl1_dir).mkdir(exist_ok=True)

            cs2_dl2_dir = os.path.join(cs2_dir,'dl2')
            Path(cs2_dl2_dir).mkdir(exist_ok=True)

            cs2_dl3_dir = os.path.join(cs2_dir,'dl3')
            Path(cs2_dl3_dir).mkdir(exist_ok=True)

            cs2_dqual_dir = os.path.join(cs2_dir,'distributions')
            Path(cs2_dqual_dir).mkdir(exist_ok=True)

            cs2_dl3_dir_eff = os.path.join(cs2_dir,'dl3_eff')
            Path(cs2_dl3_dir_eff).mkdir(exist_ok=True)

            cs2_dqual_dir_eff = os.path.join(cs2_dir,'distributions_eff')
            Path(cs2_dqual_dir_eff).mkdir(exist_ok=True)

        ## STEREO
        if stereo:
            stereo_dir = os.path.join(target_dir,'stereo')
            Path(stereo_dir).mkdir(exist_ok=True)

            stereo_dl1_dir = os.path.join(stereo_dir,'dl1')
            Path(stereo_dl1_dir).mkdir(exist_ok=True)

            stereo_dl2_dir = os.path.join(stereo_dir,'dl2')
            Path(stereo_dl2_dir).mkdir(exist_ok=True)

            stereo_dl3_dir = os.path.join(stereo_dir,'dl3')
            Path(stereo_dl3_dir).mkdir(exist_ok=True)

            stereo_dqual_dir = os.path.join(stereo_dir,'distributions')
            Path(stereo_dqual_dir).mkdir(exist_ok=True)

            stereo_dl3_dir_eff = os.path.join(stereo_dir,'dl3_eff')
            Path(stereo_dl3_dir_eff).mkdir(exist_ok=True)

            stereo_dqual_dir_eff = os.path.join(stereo_dir,'distributions_eff')
            Path(stereo_dqual_dir_eff).mkdir(exist_ok=True)
        ## R0 -> DL1
        # cs1
        aargs = iargs()
        
        aargs.config_file = config_file
        if cs1 and (target in list(dict_list_t1.keys()) ) and run_dl1:
            aargs.out_dir     = cs1_dl1_dir
            

            pool = mp.Pool(n_proc)
            pool.map(r0_dl1_1file,
                                args_maker(aargs,
                                           dict_list_t1[target]))
            pool.close()

        # cs2
        if cs2 and (target in list(dict_list_t2.keys())) and run_dl2:
            aargs.out_dir     = cs2_dl1_dir

            pool = mp.Pool(n_proc)
            pool.map(r0_dl1_1file,
                                args_maker(aargs,
                                           dict_list_t2[target]))
            pool.close()

        ## DL1 -> DL2
        aargs.models_dir  = mono_model_dir
        #cs1
        if cs1:
            if run_dl2:
                aargs.out_dir     = cs1_dl2_dir

                pool = mp.Pool(n_proc)
                pool.map(dl1_dl2_1file_1tel,
                                    args_maker(aargs,
                                               glob.glob(cs1_dl1_dir+'/*.h5')))
                pool.close()

            ## DL3
            if run_dl3:
                try:
                    aargs.input_dir   = cs1_dl2_dir
                    aargs.out_dir     = cs1_dl3_dir
                    aargs.target_name = target
                    aargs.target_ra = target_pos.ra.to_value(u.deg)
                    aargs.target_dec = target_pos.dec.to_value(u.deg)
                    aargs.irf_dir = daily_config["irf_dir"]
                    aargs.gammaness_cut_dir = None
                    dl2_dl3_1dir(aargs)
                except Exception:
                    logging.error("CS1 : DL2 > DL3 failed")

            ## efficiency cuts
            if run_dl3_eff:
                try:
                    aargs.input_dir   = cs1_dl2_dir
                    aargs.out_dir     = cs1_dl3_dir_eff
                    aargs.target_name = target
                    aargs.target_ra = target_pos.ra.to_value(u.deg)
                    aargs.target_dec = target_pos.dec.to_value(u.deg)
                    aargs.irf_dir = daily_config["irf_dir"]
                    aargs.gammaness_cut_dir = eff_cut_dir_mono
                    dl2_dl3_1dir(aargs)
                except Exception:
                    logging.error("CS1 : DL2 > DL3_eff failed")

        # cs2
        if cs2:
            if run_dl2:
                aargs.out_dir     = cs2_dl2_dir

                pool = mp.Pool(n_proc)
                pool.map(dl1_dl2_1file_1tel,
                                    args_maker(aargs,
                                               glob.glob(cs2_dl1_dir+'/*.h5')))
                pool.close()
            ## DL3
            if run_dl3:
                try:
                    aargs.input_dir   = cs2_dl2_dir
                    aargs.out_dir     = cs2_dl3_dir
                    aargs.target_name = target
                    aargs.target_ra = target_pos.ra.to_value(u.deg)
                    aargs.target_dec = target_pos.dec.to_value(u.deg)
                    aargs.irf_dir = daily_config["irf_dir"]
                    aargs.gammaness_cut_dir = None
                    dl2_dl3_1dir(aargs)
                except Exception:
                    logging.error("CS2 : DL2 > DL3 failed")

            ## efficiency cuts
            if run_dl3_eff:
                try:
                    aargs.input_dir   = cs2_dl2_dir
                    aargs.out_dir     = cs2_dl3_dir_eff
                    aargs.target_name = target
                    aargs.target_ra = target_pos.ra.to_value(u.deg)
                    aargs.target_dec = target_pos.dec.to_value(u.deg)
                    aargs.irf_dir = daily_config["irf_dir"]
                    aargs.gammaness_cut_dir = eff_cut_dir_mono
                    dl2_dl3_1dir(aargs)
                except Exception:
                    logging.error("CS2 : DL2 > DL3_eff failed")

        ## STEREO
        if stereo:
            if run_dl1:
                aargs.out_dir     = stereo_dl1_dir
                aargs.tel2_dir    = cs2_dl1_dir
                pool = mp.Pool(n_proc)
                pool.map(dl1_dl1_1file, args_maker(aargs,
                                                   glob.glob(cs1_dl1_dir+'/*.h5')))
                pool.close()       
            if run_dl2:
                aargs.out_dir   = stereo_dl2_dir
                aargs.models_dir = stereo_model_dir

                pool = mp.Pool(n_proc)
                pool.map(dl1_dl2_1file, args_maker(aargs,
                                                   glob.glob(stereo_dl1_dir+'/*.h5')))
                pool.close()
            ## DL3
            if run_dl3:
                try : 
                    aargs.input_dir   = stereo_dl2_dir
                    aargs.out_dir     = stereo_dl3_dir
                    aargs.target_name = target
                    aargs.target_ra   = target_pos.ra.to_value(u.deg)
                    aargs.target_dec  = target_pos.dec.to_value(u.deg)
                    aargs.irf_dir     = daily_config["irf_dir"]
                    aargs.gammaness_cut_dir = None
                    dl2_dl3_1dir(aargs)
                except Exception:
                    logging.error("STEREO : DL2 > DL3 failed")

            ## efficiency cuts
            if run_dl3_eff:
                try : 
                    aargs.input_dir   = stereo_dl2_dir
                    aargs.out_dir     = stereo_dl3_dir_eff
                    aargs.target_name = target
                    aargs.target_ra   = target_pos.ra.to_value(u.deg)
                    aargs.target_dec  = target_pos.dec.to_value(u.deg)
                    aargs.irf_dir     = daily_config["irf_dir"]
                    aargs.gammaness_cut_dir = eff_cut_dir_stereo
                    dl2_dl3_1dir(aargs)
                except Exception:
                    logging.error("STEREO : DL2 > DL3 failed")


        ## EXTRACT RATE DISTRIBUTION
        aargs.date_str = "{:04d}{:02d}{:02d}".format(year,month,day)
        if cs1 & run_data_qual:
            aargs.dl1_dir = cs1_dl1_dir
            aargs.dl3_dir = cs1_dl3_dir
            aargs.out_dir = cs1_dqual_dir
            extract_dl1_distributions(aargs)
        if cs2 & run_data_qual:
            aargs.dl1_dir = cs2_dl1_dir
            aargs.dl3_dir = cs2_dl3_dir
            aargs.out_dir = cs2_dqual_dir
            extract_dl1_distributions(aargs)
        if stereo & run_data_qual:
            aargs.dl1_dir = stereo_dl1_dir
            aargs.dl3_dir = stereo_dl3_dir
            aargs.out_dir = stereo_dqual_dir
            extract_dl1_distributions(aargs)       

        ## EXTRACT RATE DISTRIBUTION EFF
        aargs.date_str = "{:04d}{:02d}{:02d}".format(year,month,day)
        if cs1 & run_data_qual_eff:
            aargs.dl1_dir = cs1_dl1_dir
            aargs.dl3_dir = cs1_dl3_dir_eff
            aargs.out_dir = cs1_dqual_dir_eff
            extract_dl1_distributions(aargs)
        if cs2 & run_data_qual_eff:
            aargs.dl1_dir = cs2_dl1_dir
            aargs.dl3_dir = cs2_dl3_dir_eff
            aargs.out_dir = cs2_dqual_dir_eff
            extract_dl1_distributions(aargs)
        if stereo & run_data_qual_eff:
            aargs.dl1_dir = stereo_dl1_dir
            aargs.dl3_dir = stereo_dl3_dir_eff
            aargs.out_dir = stereo_dqual_dir_eff
            extract_dl1_distributions(aargs)  
    

    logging.info("Daily analysis ({:04d}/{:02d}/{:02d}) ended ".format(year,month,day))

if __name__ == '__main__':

    #os.environ['GAMMAPY_DATA'] = "/data/work/analysis/Daily_analysis/daily_data/"

    if len(sys.argv)>1:
        daily_config_file = sys.argv[1]
    else:
        daily_config_file = DEFAULT_CONFIG_FILE

    with open(daily_config_file) as json_data_file:
        daily_config = json.load(json_data_file)

    run_daily_ana(daily_config)




