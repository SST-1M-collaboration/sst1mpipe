#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created Feb 29 2024
"""

import glob

import astropy
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from ctapipe.io import read_table
from scipy import interpolate

from sst1mpipe.io.sst1m_event_source import SST1MEventSource

## spline aprox NSB : VAR[ADC]->SHIFT[ADC]



TCK1 = (np.array([ 0.        ,  0.        ,  0.        , 46.96969697, 60.90909091,
                   67.57575758, 74.24242424, 74.24242424, 74.24242424]),
        np.array([ -0.4367913 ,   9.08788226,  39.66927065,  54.18389959,
                   69.24518367, 104.29418015,   0.        ,   0.        ,
                   0.        ]),
        2)


TCK2=(np.array([0.        ,   0.        ,  93.74735497, 126.46232902,
                181.54405469, 256.80398702, 266.62290877, 292.98812266,
                294.71402901, 295.67072519, 298.04051867, 339.91032956,
                339.91032956]),
      np.array([ -1.74764142,  48.41280122,  68.36843193, 104.62272377,
                168.14318275, 179.34289593, 212.07400565, 214.10474285,
                215.17017626, 217.7828114 , 268.19128994,   0.        ,
                  0.        ]),
      1)

BINWIDTH = 4e-3 ## micro second

#we assume drops VS BLS (Baseline shift) to be linear :
# these numbers come from fig. 11.1 in "SiPM behaviour under continuous light"
# https://arxiv.org/abs/1910.00348
Vdrop_lin_aprox = dict({'gain_tel21': -0.00487,
                        'gain_tel22': -0.00118,
                        'PDE_tel21' : -0.00317,
                        'PDE_tel22' : -0.000663,
                        'xt_tel21'  : -0.00768,
                        'xt_tel22'  : -0.00207 })
#Sipm gain w/o Vdrop
## total pulse intergal 
Gain_0 = dict({'tel21': 21.0,
               'tel22': 23.5})


def VAR_to_shift(baseline_VAR,ntel):
    if ntel==21:
        shift = interpolate.splev(baseline_VAR, TCK1, der=0)
    elif ntel==22:
        shift = interpolate.splev(baseline_VAR, TCK2, der=0)
    else:
        print("ERROR ntel should be 21 or 22")
        shift=0
    return shift

def NSB_to_BLS(NSB_MHz,ntel):
    A = NSB_MHz*Gain_0['tel{}'.format(ntel)]*BINWIDTH
    return A/(1-Vdrop_lin_aprox['gain_tel{}'.format(ntel)]*A)

def make_drop_func(param_name,ntel):
    
    def drop_func(B_shift):
        slope   = Vdrop_lin_aprox[param_name+'_tel{}'.format(ntel)]
        #return max(slope*B_shift+1, 0)
        return slope*B_shift+1
    return drop_func

def VAR_to_Idrop(baseline_VAR,ntel):
    ## X-talk not included yet
    ## Usage :
    ## I_corr = I / I_drop
    try:
        slope   = Vdrop_lin_aprox['gain_tel{}'.format(ntel)]+ \
                  Vdrop_lin_aprox['PDE_tel{}'.format(ntel)]
    except Exception:
        print("ERROR ntel should be 21 or 22")
    return VAR_to_shift(baseline_VAR,ntel)*slope+1

def BLS_to_Idrop(baseline_shift,ntel):
    ## X-talk not included yet
    ## Usage :
    ## I_corr = I / I_drop
    try:
        slope   = Vdrop_lin_aprox['gain_tel{}'.format(ntel)]+ \
                  Vdrop_lin_aprox['PDE_tel{}'.format(ntel)]
    except Exception:
        print("ERROR ntel should be 21 or 22")
    return baseline_shift*slope+1


def VAR_to_Gdrop(baseline_VAR,ntel):
    ## Usage :
    ## G_corr = G / G_drop
    try:
        slope   = Vdrop_lin_aprox['gain_tel{}'.format(ntel)]
    except Exception:
        print("ERROR ntel should be 21 or 22")

    return VAR_to_shift(baseline_VAR,ntel)*slope+1


def get_simple_nsbrate(bs_shift,gain,binlenght=4e-3):
    rate = bs_shift/gain / binlenght
    return rate

def gain_drop_th(nsb_rate, cell_capacitance=85. * u.fF, bias_resistance=2.4 * u.kohm):
    return 1 - 1 / (1 + (nsb_rate * cell_capacitance * bias_resistance).to_value(1))


    
def BLS_to_NSB(baseline_shift,ntel,gain=None):
    if gain is None:
        gain = Gain_0['tel{}'.format(ntel)]
    Gdrop = baseline_shift * Vdrop_lin_aprox['gain_tel{}'.format(ntel)]+1
    NSB_rate = get_simple_nsbrate(baseline_shift,
                                  gain * Gdrop)
    return NSB_rate

def BLS_to_NSB_photon_rate(baseline_shift,ntel,gain=None):
    PDE_0 = 0.35 ## aprox. TODO : get some better estimation
    if gain is None:
        gain = Gain_0['tel{}'.format(ntel)]

    slope   = Vdrop_lin_aprox['gain_tel{}'.format(ntel)]+ \
              Vdrop_lin_aprox['PDE_tel{}'.format(ntel)]

    Idrop = baseline_shift*slope+1
    NSB_rate = get_simple_nsbrate(baseline_shift,
                                  gain * Idrop) / PDE_0

    return NSB_rate

def VAR_to_NSB(baseline_VAR,ntel,gain=None):
    if gain is None:
        gain = Gain_0['tel{}'.format(ntel)]
    NSB_rate = BLS_to_NSB(VAR_to_shift(baseline_VAR,ntel),
                          ntel,
                          gain=gain)
    return NSB_rate

def VAR_to_NSB_photon_rate(baseline_VAR,ntel,gain=None):
    if gain is None:
        gain = Gain_0['tel{}'.format(ntel)]
    NSB_rate = BLS_to_NSB_photon_rate(VAR_to_shift(baseline_VAR,ntel),
                                      ntel,
                                      gain=gain)
    return NSB_rate

def find_dark_files(data_dir):
    dark_files = []
    dark_names=['dark','DARK','drak','DRAK',"Dark","Drak"]
    file_list=glob.glob(data_dir+"/*.fits.fz")

    for ii,file_path in enumerate(file_list):
        #run_number = int(file_path.split("_")[-1].split('.')[0])
        try :
            f = astropy.io.fits.open(file_path)
            target = f[2].header['TARGET'] 
        except Exception:
            print('failed reading {}'.format(file_path))
            continue
        if target in dark_names:
            dark_files.append(file_path)
        f.close()
    return dark_files


def get_dark_baseline(filename,max_evt=500,event_type=8):
        
        raw_baselines  = [ [] for ii in range(1296)]
       
        data_stream = SST1MEventSource([filename],
                                       max_events=max_evt)
        for ii,event in enumerate(data_stream):
                    tel = event.sst1m.r0.tels_with_data[0]
                    r0data = event.sst1m.r0.tel[tel]

                    if r0data._camera_event_type.value==event_type:
                        for pix in range(1296):
                            raw_baselines[pix].append(r0data.adc_samples[pix,:50])
        raw_baselines  = np.array(raw_baselines)
        
        return raw_baselines.mean(axis=(1,2))

############
 #### tools to plot data :
############

def get_ped_table(file_list):
    
    bline_table = None
    
    for dl1file in sorted(file_list):
        if bline_table is None :
            try:
                bline_table = read_table(dl1file,
                                         '/dl1/monitoring/telescope/pedestal')
            except Exception:
                print("pedestal not found in {}".format(dl1file))
        else :
            try:
                bline_table = astropy.table.vstack([bline_table,
                                                   read_table(dl1file,
                                                              '/dl1/monitoring/telescope/pedestal')])
            except Exception:
                print("pedestal not found in {}".format(dl1file))
    return bline_table

def get_ped_table_low_res(file_list):
    
    bline_table = None
    
    for dl1file in sorted(file_list):
        if bline_table is None :
            try:
                bline_table = read_table(dl1file,
                                         '/dl1/monitoring/telescope/pedestal')[-1]
            except Exception:
                print("pedestal not found in {}".format(dl1file))
        else :
            try:
                bline_table = astropy.table.vstack([bline_table,
                                                   read_table(dl1file,
                                                              '/dl1/monitoring/telescope/pedestal')[-1]])
            except Exception:
                print("pedestal not found in {}".format(dl1file))
    return bline_table


def plot_average_nsb_VS_time(ped_table,ntel,ax=None, color='blue',label='label'):
    NSB = VAR_to_NSB(ped_table['pedestal_charge_std'].mean(axis=1)**2,ntel)
    Dates = [Time(t,scale='utc',format='unix').to_datetime() for t in  ped_table['pedestal_sample_time']]
    if ax is None:
        f,ax = plt.subplots(figsize=(10,5))
    ax.plot(Dates,NSB,'.',label=label, color=color)
    plt.xlabel('Time')
    plt.ylabel('NSB [MHz]')
    return ax

def plot_average_nsb_photon_rate_VS_time(ped_table,ntel,ax=None):
    NSB = VAR_to_NSB_photon_rate(ped_table['pedestal_charge_std'].mean(axis=1)**2,ntel)
    Dates = [Time(t,scale='utc',format='unix').to_datetime() for t in  ped_table['pedestal_sample_time']]
    if ax is None:
        f,ax = plt.subplots(figsize=(10,5))
    ax.plot(Dates,NSB,'.',label='tel {}'.format(ntel%20))
    plt.xlabel('Time')
    plt.ylabel('NSB photons [MHz]')
    return ax
