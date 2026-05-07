#!/usr/bin/env python3
"""
Created June 01 2024
"""

import glob

import astropy.units as u
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.time import Time
from ctapipe.io import read_table
from scipy.optimize import curve_fit

from sst1mpipe.io import load_config, load_dl1_sst1m
from sst1mpipe.performance import get_mc_info, get_weights
from sst1mpipe.performance.spectra import DAMPE_P_He_SPECTRUM

DEFAULT_CUTS = {}
DEFAULT_CUTS["zenith"]             = [0,70]
DEFAULT_CUTS["livetime"]           = [0,5]
DEFAULT_CUTS["ped_fraction"]       = [0,0.05]
DEFAULT_CUTS["tc_raised_fraction"] = [0,1]
DEFAULT_CUTS["NSB"]                = [0,600]
DEFAULT_CUTS["MC_rate_ratio"]      = [0.7,1.3]
DEFAULT_CUTS["failed_fit"]         = [0,0.5]


def get_slow_data_table(dl3_file,file_radical='DIGICAM',root_dir='/net'):
    dl3 = fits.open(dl3_file)
    
    datestr = str(dl3[1].header["OBS_ID"])[:8]
    oyear   = str(dl3[1].header["OBS_ID"])[:4]
    omonth  = str(dl3[1].header["OBS_ID"])[4:6]
    oday    = str(dl3[1].header["OBS_ID"])[6:8]
    try:
        itel = int(dl3[1].header['TELESCOP'][-2:])%20
    except (IndexError, KeyError, ValueError, TypeError):
        print("This dont work for setereo DL3")
        itel = None
    
    tstart = Time(dl3[1].header['TSTART'],format='unix',scale='utc')
    tstop  = Time(dl3[1].header['TSTOP'],format='unix',scale='utc')
    
    files = glob.glob(f'{root_dir}/cs{itel}/data/aux/{oyear}/{omonth}/{oday}/SST1M_TEL{itel}/{file_radical}{itel}_{datestr}_*.fits')
    table = vstack([ Table(fits.open(f)[1].data) for f in files])
    t_mask = (table['TIMESTAMP']>tstart.unix*1000) & (table['TIMESTAMP']<tstop.unix*1000)
    return table[t_mask]

def get_bad_intervals(table,gapmax=2000,tpas=60000):
    sat_mask = table['err_udpktarena_nfull']<gapmax
    if (~sat_mask).sum()==0:
        return []
    tt = t0  = table['TIMESTAMP'][~sat_mask][0]
    NGTI =[]

    for _,t in enumerate(table['TIMESTAMP'][~sat_mask]):
        if np.abs(t-tt) > tpas:
            NGTI.append([t0,tt])
            t0=t
        tt=t
    NGTI.append([t0,t])
    return NGTI



def plot_rates_from_slow_data(dl3_files):
    _,ax =plt.subplots(figsize=(15,6))

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())

    for dl3_file in sorted(dl3_files):
        obs_id = dl3_file.split("_")[-2]
        table = get_slow_data_table(dl3_file)
        table_dsc = get_slow_data_table(dl3_file,file_radical='DigicamSlowControl')

        ax.plot([Time(tt/1000.,format='unix').to_datetime() for tt in table['TIMESTAMP'] ] ,
                     table['rate_swat_arrevent_nreads'],
                     '.',
                     alpha=.3,
                     label=f'obs_id : {obs_id}',
                     #label='rate_swat_arrevent_nreads (rate of array events successfully read from SWAT)'
                     )


        ax.plot([Time(tt/1000.,format='unix').to_datetime() for tt in table_dsc['TIMESTAMP'] ] ,
                     table_dsc['triggerStatus'][:,3],
                     '-',
                     alpha=.3,
                     color='black',
                     )#label='rate from Digicam slow control')



        NGTI = get_bad_intervals(table)
        for interval in NGTI:
            if len(NGTI)>0:
                ax.axvspan(Time(interval[0]/1000.,format='unix').to_datetime(),
                           Time(interval[1]/1000.,format='unix').to_datetime(),
                           alpha=0.2, color='red')
    ax.set_yscale('log')
    ax.set_ylabel('Rate [Hz]')
    plt.legend()
    return ax

def get_MC_dist_mono(tel_setup,
                     MC_proton_file,
                     config_file,
                     zenith,
                     Q_bins):
    
    mc_config = load_config(config_file)
    dl2_mc_proton = load_dl1_sst1m(MC_proton_file,
                                  tel=tel_setup,
                                  config=mc_config,
                                  table='astropy')

    mc_info_proton = get_mc_info(MC_proton_file, config=mc_config)

    dl2_mc_proton  = get_weights(dl2_mc_proton, 
                                 mc_info=mc_info_proton, 
                                 obs_time=1*u.s, 
                                 target_spectrum=DAMPE_P_He_SPECTRUM)
    
    rates_mc,_ = np.histogram(dl2_mc_proton["camera_frame_hillas_intensity"],
                                       weights=dl2_mc_proton["weight"],
                                       bins=Q_bins)
    diff_rates_mc = rates_mc/(Q_bins[1:]-Q_bins[:-1])/np.cos(zenith*u.deg)


    return diff_rates_mc

def write_MC_dist(tel_setup,
                  MC_proton_file,
                  config_file,
                  zenith,
                  outdir,
                  Q_bins = np.geomspace(10,1e5,300)):
    
    res_dict = {}
    res_dict["low"] = Q_bins[:-1]
    res_dict["high"] = Q_bins[1:]
    res_dict["center"] = (Q_bins[:-1] + Q_bins[1:]) / 2
        
    if tel_setup=='stereo':
        tel_setups = ['tel_001','tel_002']
    else:
        tel_setups = [tel_setup]
    for tel in tel_setups:
        diff_rates_mc = get_MC_dist_mono(tel,
                                         MC_proton_file,
                                         config_file,
                                         zenith,
                                         Q_bins)

        res_dict["diff_rate_zenCorected_"+tel] = diff_rates_mc
        
    res_df = pd.DataFrame(res_dict)
    outfile = outdir+f'/MC_{tel_setup}_intensity_hist.h5'
    res_df.to_hdf(outfile,'intensity_hist')

def getmask(key,sel_dict,DQ_table):
    mask = (DQ_table[key]>=sel_dict[key][0]) & (DQ_table[key]<=sel_dict[key][1])
    return np.array(mask)

def make_selection(DQ_table,
                   sel_dict=DEFAULT_CUTS):
    
    flag_array = np.ones(DQ_table.shape[0],dtype=bool)
    for key in sel_dict.keys():
        flag_array = flag_array & getmask(key,sel_dict,DQ_table)
    return np.array(flag_array)


def make_DQ_table(tel_setup,
                  file_list,
                  mc_hist_file,
                  outdir,
                  Q_min=400,
                  Q_max=4000,
                  sel_dict=DEFAULT_CUTS):
    def lin(x,a,b):
        return a*x+b
    if tel_setup=='stereo':
        tel_setups = ['tel_021','tel_022']
    else:
        tel_setups = [tel_setup]
    for tel in tel_setups:
        res_dict = {}
        res_dict["obs_id"]             = []
        res_dict["zenith"]             = []
        res_dict["livetime"]           = []
        res_dict["ped_fraction"]       = []
        res_dict["tc_raised_fraction"] = []
        res_dict["NSB"]                = []
        res_dict["MC_rate_ratio"]      = []
        res_dict["failed_fit"]         = []
        res_dict["qual_flag"]          = []

        MC_hist  = pd.read_hdf(mc_hist_file)
        rates_mc = MC_hist[f'diff_rate_zenCorected_tel_00{tel[-1]}']
        Q_bins_c = MC_hist['center']

        fitmask = (Q_bins_c>Q_min) & (Q_bins_c<Q_max)
        mc_fit, mc_err = curve_fit(lin,
                                   np.log10(Q_bins_c[fitmask]),
                                   np.log10(rates_mc[fitmask])
                                   )

        def shifted_lin(x,shift):
            return lin(x-shift,*mc_fit)


        for _ ,filename in enumerate(file_list):
            if tel in filename:
                int_hist = read_table(filename,'intensity_hist')
                zenith   = read_table(filename,'zenith')[0][0]
                obs_id   = filename.split('_')[-1].split('.')[0]
                res_dict["obs_id"].append(int(obs_id))
                res_dict["zenith"].append(zenith)
                res_dict["livetime"].append(read_table(filename,'t_elapsed')[0][0])
                res_dict["ped_fraction"].append(read_table(filename,'survived_pedestal_frac')[0][0])
                res_dict["tc_raised_fraction"].append(read_table(filename,'recleaned_fraction')[0][0])
                res_dict["NSB"].append(read_table(filename,'NSB')[0][0])
                res_dict["qual_flag"].append(0)
                try:
                    Q_bins = np.append(int_hist['low'],int_hist['high'][-1])
                    Q_bins_c = (Q_bins[1:]+Q_bins[:-1])/2.
                    fitmask = (Q_bins_c>Q_min) & (Q_bins_c<Q_max)
                    no_z_mask = fitmask & (int_hist['diff_rate']>0)
                    popt,corr= curve_fit(shifted_lin,
                                         np.log10(int_hist['center'][no_z_mask]),
                                         np.log10(int_hist['diff_rate'][no_z_mask]),
                                         )  
                    res_dict["MC_rate_ratio"].append(10**popt[0]/np.cos(zenith*u.deg).to_value())
                    res_dict["failed_fit"].append(0)
                except (RuntimeError, ValueError, TypeError):
                    res_dict["MC_rate_ratio"].append(0)
                    res_dict["failed_fit"].append(1)       
                    print("fit failed",filename)

        res = pd.DataFrame(res_dict)
        res["qual_flag"] = make_selection(res,sel_dict=sel_dict)
        outfile = outdir+f'/DQ_table_{tel_setup}.h5'
        res.to_hdf(outfile,f'DQ_table_{tel}')
    
    return

