#!/usr/bin/env python3
"""
Created on Wed Mar 23 16:36:22 2022

@author: zrdz
"""


import argparse
import datetime
import os

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

#from spe_mc_v6 import mes_fitter_mc
from astropy.io import ascii as asc

#from cts_core.camera import Camera
#from digicampipe.instrument import geometry
#from digicampipe.io.event_stream import event_stream, add_slow_data
from ctapipe.io import EventSource
from PulseTemplate import GetTemplate3
from scipy.interpolate import UnivariateSpline
from scipy.optimize import minimize

from sst1mpipe.calib import get_default_calibration
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.utils import get_cam_geom


class shape_maker:
    def __init__(self,
                 day            = 26,
                 month          = 3,
                 year           = 2022,
                 data_path      = "/net/",
                 tel            = 1,
                 first_file_n   = 1,
                 n_files        = 1,
                 plot_dir       = './',
                 bshift         = 0,
                 sample_size    = 4*u.ns,
                 max_evt        = None,
                 dark_baselines = None,
                 pix            = 321,
                 MC_filename    = None):
        if MC_filename is None:
            self.isMC = False
        else:
            self.isMC = True

        self.tel=tel
        self.pix = pix
        date_str = f'{year:04d}{month:02d}{day:02d}'
        self.sample_size = sample_size

        self.data_path = data_path
        self.files_path   = os.path.join(data_path,
                                         f'cs{tel}',
                                         'data',
                                         'raw',
                                         f'{year:04d}',
                                         f'{month:02d}',
                                         f'{day:02d}',
                                         f'SST1M{tel}')
        self.first_file_n = first_file_n
        self.n_files      = n_files
        self.plot_dir     = plot_dir
        self.max_evt      = max_evt
        self.bshift       = bshift

        self.dark_baselines = dark_baselines


        self.geom = get_cam_geom(20+self.tel)

        self.n_pixels = self.geom.n_pixels
        self.pixels   = np.arange(self.n_pixels)

        self.file_list = []
        self.Q7ns_array=np.array([])
        self.Qtot_array=np.array([])
        if self.isMC:
            self.file_list.append(MC_filename)
        else:
            for n in range(first_file_n, first_file_n+n_files):
                filepath = os.path.join(self.files_path,f'SST1M{self.tel}_{date_str}_{n:04d}.fits.fz')
                if os.path.isfile(filepath):
                    self.file_list.append(filepath)
                else:
                    print(f'file {filepath} not found')
            if len(self.file_list) == 0:
                print("Warning : no files")

        self.res = None

    def make_shape(self,n_rebin=4,ManualShift = 0):

        def get_barry_tom(wf):
                tmax = np.argmax([ wf[ii:ii+3].sum() for ii in range(len(wf)-3) ])
                return (wf[tmax:tmax+3]*np.arange(tmax,tmax+3)/wf[tmax:tmax+3].sum()).sum()
        def rebin(samples,n_rebin):
            rb_samples=[]
            samples = samples.astype(float)
            for sub_sample in range((len(samples)*n_rebin)-n_rebin):
                rb_samples.append(samples[sub_sample//n_rebin]+(sub_sample%n_rebin)/n_rebin*(samples[(sub_sample//n_rebin)+1]-samples[(sub_sample//n_rebin)]))
                #print(type((samples[sub_sample//n_rebin]-samples[(sub_sample//n_rebin)+1])))
            for sub_sample in range(n_rebin):
                rb_samples.append(samples[len(samples)-1])
            return rb_samples

        def likelyhood_fct_spe(fct,wf,shift):
            like = 0
            xx = np.arange(50)

            like = np.nansum(fct(xx+shift)*wf)
            return -like


        def splitsum(t1,t2,shift):
            tmpf = UnivariateSpline(np.array(range(len(t2)))+shift,t2,ext=2,k=1,s=0)
            ### DAm bug in UnivariateSpline.intergal assume ext=0 even if not..

            #~ ttmp = [tmpf.integral(t,t+1) + t1_val for t,t1_val in enumerate(t1)]
            ttmp =[]
            pp =0
            # pp0 =0
            for t,t1_val in enumerate(t1):
                #print(tmpf.integral(t,t+1))
                try :
                    pp = tmpf(t+1)
                    # pp0 = tmpf(t)
                    tttmp = tmpf.integral(t,t+1) + t1_val
                except ValueError:
                    tttmp = t1_val + pp

                ttmp.append(tttmp)
            #plt.plot(tmpf(np.arange(0,60*4+50,1)))
            #plt.show()
            return ttmp

        pixs_interp_norm,T0 = GetTemplate3()

        pulseshape = np.zeros(50*n_rebin)
        maxref_0 = 14.2
        shifts =[]
        LLs=[]
        if self.isMC:
            data_stream = EventSource(
                self.file_list[0],
                max_events=self.max_evt
                )
            mask_low_el_noise = np.ones(self.n_pixels,dtype=bool)
        else:
            data_stream = SST1MEventSource(
                filelist=self.file_list,
                disable_bar = True,
                max_events=self.max_evt
                )
            calib_param = get_default_calibration(self.tel+20)[0]
            mask_low_el_noise = np.array(calib_param['sigma_el']<5)
            mask_low_el_noise = mask_low_el_noise & np.array(calib_param['calib_flag']==1)

        tot_evts = 0

        for ii,event in enumerate(data_stream):
            #for tel in event.sst1m.r0.tels_with_data:
                tel = self.tel+20
                if self.isMC:
                    r0data = event.r0.tel[tel]
                else:
                    r0data = event.sst1m.r0.tel[tel]

                if ii==0 and not self.isMC:
                    self.T0 = r0data.local_camera_clock/1e9
                    self.start_date = datetime.datetime.fromtimestamp(self.T0)


                # sel_wfs = []




                if self.isMC:
                    mcdata = event.mc.tel[tel]
                    Qsum = (r0data.adc_samples[0].T[20:35] - mcdata.pedestal/50.+self.bshift).sum(axis=0)
                    wfs = (r0data.adc_samples[0].T - mcdata.pedestal/50.+self.bshift).T

                else:
                    if r0data._camera_event_type.value==8:
                        if self.dark_baselines is None:
                            Qsum = (r0data.adc_samples.T[20:35]+self.bshift - r0data.digicam_baseline).sum(axis=0)
                            wfs = (r0data.adc_samples.T+self.bshift - r0data.digicam_baseline).T

                        else:
                            Qsum = (r0data.adc_samples.T[20:35] - self.dark_baselines ).sum(axis=0)
                            wfs = (r0data.adc_samples.T - self.dark_baselines).T
                    else:
                        continue


                if self.pix is None:
                    if self.tel==1:
                        mask1pe = (Qsum>13)*(Qsum<24)
                    else:
                        mask1pe = (Qsum>16)*(Qsum<27)
                    for pix in self.pixels[mask1pe & mask_low_el_noise]:
                        wf = wfs[pix]



                        btom = get_barry_tom(wf)
                        if btom<20 or btom>35:
                            continue

                        QQ = wf[round(btom)-3:round(btom)+4].sum()
                        QQtot = wf[round(btom)-6:round(btom)+14].sum()
                        if QQ < 13:
                            continue
                        if wf[:20].sum(axis=0) >14:
                            continue

                        if wf[35:].sum(axis=0) >14:
                            continue

                        try:
                            P0=T0-btom
                            res1  = minimize( lambda x: likelyhood_fct_spe(pixs_interp_norm,wf,x)   ,[P0] , method='BFGS', tol=1e-12)
                            tom0  = T0-res1.x[0]
                            # LL = 2*(np.log(res1.fun) - np.log(likelyhood_fct_spe(pixs_interp_norm,wf,P0)))
                            # print(tom0,'diff :',res1.x[0]-P0)
                            LLs.append(res1.x[0]-P0)
                            shift = (maxref_0*n_rebin)-tom0*n_rebin

                            xx = np.linspace(0,50,200)



                        except Exception:
                            print("fit failed : skippyng")
                            continue

                        if np.random.random() >0.9 and False:
                        #if QQ<10:
                            f,ax = plt.subplots()
                            ax.plot(wf)
                            ax.plot(xx,pixs_interp_norm(xx),alpha=.3)
                            ax.plot(xx,pixs_interp_norm(xx+T0-btom),'--')
                            ax.plot(xx,pixs_interp_norm(xx+T0-tom0))
                            #plt.show()

                        shape = rebin((wf),n_rebin)
                        tot_evts +=1
                        shifts.append(shift)
                        pulseshape = splitsum(pulseshape,shape,shift)
                        self.Q7ns_array = np.append(self.Q7ns_array,QQ)
                        self.Qtot_array = np.append(self.Qtot_array,QQtot)
                else:
                    if Qsum[self.pix]>15 and Qsum[self.pix]<30:
                        wf = wfs[self.pix]
                        btom = get_barry_tom(wf)
                        if btom<20 or btom>35:
                            continue

                        shift = (maxref_0*n_rebin)-btom*n_rebin
                        QQ = wf[round(btom)-3:round(btom)+4].sum()
                        shape = rebin((wf/QQ),n_rebin)

                        tot_evts +=1
                        shifts.append(btom-maxref_0)
                        pulseshape = splitsum(pulseshape,shape,shift)
                        #print(np.argmax(ttt))
        self.tot_evts = tot_evts
        t_0 = np.argmax(pulseshape)/float(n_rebin)
        print(f'T_0 : {t_0} ns')
        rebins = np.linspace(0,50-(1/n_rebin),50*n_rebin)

        ##ManualShift is here to set the TOM at 30 nanosec (or any val)
        ## Qamp is the mean amplitude/Q ratio

        ManualShift = 0
        print("shft",ManualShift)
        print("Tstd :: ",np.std(shifts))
        print('Tot Evts :',self.tot_evts)
        # ~ pulse_fct = UnivariateSpline(rebins+ManualShift,np.array(pulseshape)*QAmp/np.max(pulseshape),ext=3,k=3,s=0)
        self.pulse_fct    = UnivariateSpline(rebins+ManualShift,np.array(pulseshape),ext=3,k=3,s=0)
        self.pulse_fct_ns = UnivariateSpline((rebins+ManualShift)*self.sample_size.to_value('ns'),
                                              np.array(pulseshape),ext=3,k=3,s=0)
        self.fit_vs_b = np.array(LLs)

        return








#########################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--plot_dir',   type=str, default = './SPE_plots_mc/' )
    args = parser.parse_args()



    shaper_t1 = shape_maker(day          = 18,
                            month        = 3,
                            year         = 2024,
                            data_path    = "/net/",
                            tel          = 1,
                            first_file_n = 738,
                            n_files      = 2,
                            max_evt      = 1000,
                            plot_dir     = args.plot_dir,
                            pix          = None)
                          # dark_baselines=bs_t1_dark.raw_baselines.mean(axis = (1,2)))
    shaper_t1.make_shape()



    shaper_t2 = shape_maker(day          = 18,
                            month        = 3,
                            year         = 2024,
                            data_path    = "/net/",
                            tel          = 2,
                            first_file_n = 706,
                            n_files      = 2,
                            max_evt      = 1000,
                            plot_dir     = args.plot_dir,
                            pix          = None)
                            # dark_baselines=bs_t2_dark.raw_baselines.mean(axis = (1,2)))
    shaper_t2.make_shape()


    shaper_mc = shape_maker(MC_filename="./simtel/biascurve_run663_TEL1_nsb2_dsum260.simtel.gz",
                          plot_dir     = args.plot_dir,
                          max_evt      = 1000,
                          pix          = 321)
    #shaper_mc.make_shape()



    pt = asc.read('./pulse_template_TEL1_3.txt')

    xx = np.linspace(0,50,400)

    f,ax = plt.subplots(figsize=(12,5))
    ax.plot(xx*4,shaper_t1.pulse_fct(xx)/shaper_t1.pulse_fct(xx).max(),label='Data Tel 1')
    ax.plot(xx*4,shaper_t2.pulse_fct(xx)/shaper_t2.pulse_fct(xx).max(),label='Data Tel 2')

    #ax.plot(xx*4,shaper_mc.pulse_fct(xx)/shaper_mc.tot_evts*gain,label='MC reco')

    ax.plot(pt['col1']+9.75*4,pt['col2'],label='MC input (template)')
    ax.legend()
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("ADC count")
    ax.set_title(' One Ph.e. Normalised Pulse shape')
    plt.grid()

    f,ax = plt.subplots(figsize=(12,5))
    ax.plot(xx*4,shaper_t1.pulse_fct(xx)/shaper_t1.tot_evts,label='Data Tel 1')
    ax.plot(xx*4,shaper_t2.pulse_fct(xx)/shaper_t2.tot_evts,label='Data Tel 2')

    #ax.plot(xx*4,shaper_mc.pulse_fct(xx)/shaper_mc.tot_evts*gain,label='MC reco')

    ax.plot(pt['col1']+9.75*4,pt['col2']*5.7,label='MC input (template)')
    ax.legend()
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("ADC count")
    ax.set_title(' One Ph.e. Pulse shape')
    plt.grid()


    plt.show()
    # ax.plot(pt['col1'],pt['col3'])






