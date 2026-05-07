#!/usr/bin/env python3
"""
Created on Wed Mar 23 16:36:22 2022
@author: TT 



    This script is ment to produce h5 table containing calibration parameters.
    This is done fitting mes spectrum using dark run data.
    it need few tens of thousan events (2 raw data file is probably enough)



    Parameters
    ----------
    plot_dir: string
        Path to a directory to save plots
    data_path: string
        Root directory of the raw data ('/net' on calculus)
    day: int
        Day of the data you want to use
    month: int
        Month of the data you want to use
    year: int
        Year of the data you want to use

    first_file: int
        number of the first file to use
    tel: int
        number of the telescope (1 or 2)
    nfile: int
        number of file to use starting from first_file 
    save_dir: string
        Path to a directory to save h5 files
    save_plot: boolean
        Save plots in  plot_dir if True
    n_proc: int
        number of core to use if run on a server with multiple proc

    Returns
    -------

    usage exemple :
    python dark_run_mes_fitter.py --year 2023 --month 2 --day 28 --tel first --2_file 129 --n_files 2 --save_plot False --save_dir ./calib

"""

import argparse
import multiprocessing as mp
import os

# from ctapipe.instrument import CameraGeometry
# from ctapipe.image import hillas_parameters, tailcuts_clean
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.special as scsp
import scipy.stats as scst
from ctapipe.instrument import SubarrayDescription
from ctapipe.visualization import CameraDisplay
from iminuit import Minuit
from pkg_resources import resource_filename

from sst1mpipe.io.sst1m_event_source import SST1MEventSource

#from iminuit.cost import LeastSquares


class mes_fitter:
    """
    class desinged to produce h5 table containing calibration parameters. 
    """
    def __init__(self,
                 day          = 5,
                 month        = 4,
                 year         = 2023,
                 data_path    = "/net/",
                 tel          = 1,
                 first_file_n = 154,
                 n_files      = 1,
                 plot_dir     = './',
                 save_dir     = './calib_h5/',
                 max_evt      = 100000,
                 peak_search_window_width = 5,
                 dark_baselines = None):
        """"
        Parameters
        ----------
        day: int
            Day of the data you want to use
        month: int
            Month of the data you want to use
        year: int
            Year of the data you want to use
        data_path: string
            Root directory of the raw data ('/net' on calculus)
        first_file: int
            number of the first file to use
        tel: int
            number of the telescope (1 or 2)
        nfile: int
            number of file to use starting from first_file 
        plot_dir: string
            Path to a directory to save plots
        save_dir: string
            Path to a directory to save h5 files
        max_evt: int
            Max number of event used
        dark_baselines : 1d array
            Array of dark baselines value for each pixels. If none digicam baseline will be used
        """

        self.tel        = tel
        self.date_str   = f'{year:04d}{month:02d}{day:02d}'
        self.save_dir   = save_dir
        self.data_path  = data_path
        self.files_path  = os.path.join(data_path,
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

        self.dark_baselines = dark_baselines


        subarray_file = resource_filename('sst1mpipe',
                                          os.path.join('data',
                                                    'sst1m_array.h5'
                                                    )
                                          )
        subarray = SubarrayDescription.from_hdf(subarray_file, focal_length_choice="EQUIVALENT")

        self.geom = subarray.tel[21].camera.geometry


        self.n_pixels = self.geom.n_pixels
        self.pixels   = np.arange(self.n_pixels)

        self.save_str = self.date_str+f"_{first_file_n:04d}_{first_file_n+n_files-1:04d}_Tel{tel}_"

        self.file_list = []
        for n in range(first_file_n, first_file_n+n_files):
            filepath = os.path.join(self.files_path,f'SST1M{tel}_{self.date_str}_{n:04d}.fits.fz')
            if os.path.isfile(filepath):
                self.file_list.append(filepath)
            else:
                print(f'file {filepath} not found')
        if len(self.file_list) == 0:
            print("Warning : no files")

        self.res = None
        self.peak_search_window_width = peak_search_window_width
        self.n_sample = 50
        self.samples=np.arange(self.n_sample)
        self.param_names = [ 'dcr','X_talk', 'gain', 'sigma_pe', 'sigma_el',"x_el"]

    ##############################################

    def set_file_list(self,file_list):
        self.file_list = file_list


    def get_histograms(self):
        """
        read raw files and and load ADC histograms
        """

        n_bins_adcmax = 51
        n_bins_adcsum = 176


        centers_adcsum = np.linspace(-50,125,n_bins_adcsum)
        centers_adcmax = np.linspace(-0,50,n_bins_adcmax)

        Qsum_hist      =  np.zeros([self.n_pixels,n_bins_adcsum])
        Qmax_hist      =  np.zeros([self.n_pixels,n_bins_adcmax])




        tot_evts = 0
        print("starting. reading data. Loading histograms.")

        data_stream = SST1MEventSource(
            self.file_list,
            disable_bar = True,
            max_events=self.max_evt
            )

        for ii,event in enumerate(data_stream):
                tel = event.sst1m.r0.tels_with_data[0]
                r0data = event.sst1m.r0.tel[tel]


                if r0data._camera_event_type.value==8:
                    tot_evts +=1

                    if self.dark_baselines is None:
                        wfs = (r0data.adc_samples.T - r0data.digicam_baseline).T

                        #sums = convolve1d(
                        #    wfs[:,:],
                        #    np.ones(self.peak_search_window_width), axis=1, mode="nearest"
                        #    )
                        #w_start = np.argmax(sums[:, self.peak_search_window_width//2:-(self.peak_search_window_width//2)],
                        #                    axis=1)

                        sums = np.array([np.convolve(wf,
                                                     np.ones(self.peak_search_window_width),
                                                     mode="full")
                                        for wf in wfs])
                        w_start = np.argmax(sums[:, self.peak_search_window_width-1:-self.peak_search_window_width+1],
                                            axis=1)



                        Qsum = np.array([wfs[ii,w_start[ii]:w_start[ii]+self.peak_search_window_width].sum() for ii in self.pixels])

                        Qmax = (r0data.adc_samples.T       - r0data.digicam_baseline).max(axis=0)
                        #Qsum_out = np.array([wfs[ii,:w_start[ii]].sum()+ \
                        #                     wfs[ii,w_start[ii]+self.peak_search_window_width:].sum() for ii in self.pixels])
                    else:
                        #Qsum = (r0data.adc_samples.T[-15:] - self.dark_baselines ).sum(axis=0)
                        #Qmax = (r0data.adc_samples.T       - self.dark_baselines ).max(axis=0)
                        pass

                    i_to_fill_adcsum = np.searchsorted(centers_adcsum[1:-1], Qsum)
                    i_to_fill_adcmax = np.searchsorted(centers_adcmax[1:-1], Qmax)

                    for pix in self.pixels:
                        Qsum_hist[pix][i_to_fill_adcsum[pix]] = Qsum_hist[pix][i_to_fill_adcsum[pix]] +1
                        Qmax_hist[pix][i_to_fill_adcmax[pix]] = Qmax_hist[pix][i_to_fill_adcmax[pix]] +1



        print(f"{ii+1} evts proceeded")
        print(f"{tot_evts} evts in histogram")
        self.tot_evts = tot_evts
        self.binwidth = centers_adcsum[1]-centers_adcsum[0]
        self.centers_adcsum = centers_adcsum
        self.Qsum_hist      = Qsum_hist
        self.centers_adcmax = centers_adcmax
        self.Qmax_hist      = Qmax_hist

        return(centers_adcsum,Qsum_hist,centers_adcmax,Qmax_hist)

    ##############################################

    def spe_spectrum_function__s(self,x, ll, xt, g, sigma_pe, sigma_el):
        """
        single phot-electron spectrum fuction x -> MES(x) 

        Parameters
        ----------
        x: integer
            ADC counts
        ll: float
            Expected value of the poisson law (avergae thermal photon number in 60ns)
        xt : float
            Cross talk probability
        g: float
            gain : Averaged number of ADC counts produced by one photo-electron
        sigma_pe: float
            st.d. of the ADC counts produced by one photo-electron
        sigma_el: float
            st.d. of the ADC counts produced by the electronic noise
        """

        STP = np.sqrt(2 * np.pi)
        x0 = -ll * g/(1-xt)

        def single_gauss(x, xn, sigma_n):
            return np.exp(-1/2 *((x-xn)/sigma_n)**2) / (sigma_n * STP)

        def sigma_n(n, sigma_pe, sigma_el):
            return np.sqrt(n * sigma_pe**2 + sigma_el**2)
        def G_poisson(n,theta,l_lambda): ## folowing notation from On the Generalized Poisson Distribution (2006) lambda = muXT!
            a1 = theta*(theta + n*l_lambda)**(n-1) / scsp.factorial(n)
            a2 = np.exp(-theta - n*l_lambda)
            return a1*a2

        def spe_spectrum_function_(x, ll, xt, g, x0, sigma_pe, sigma_el):

            S = np.zeros_like(x)
            for n in range(0, 6):
                if n==0:
                    S+= single_gauss(x, x0, sigma_el) * G_poisson(0,ll,xt)
                else:
                    xn =  n * g
                    _sigma_n = sigma_n(n, sigma_pe, sigma_el)

                    S += (
                         G_poisson(n,ll,xt)*
                         single_gauss(x, xn, _sigma_n)
                         )

            return  S

        return spe_spectrum_function_(x, ll, xt, g, x0, sigma_pe, sigma_el) * self.tot_evts* self.binwidth

    ##############################################
    ##############################################
    ##############################################

    def spe_spectrum_function_sss(self,x, ll, xt, g, sigma_pe, sigma_el,x_el):
        """
        single phot-electron spectrum fuction x -> MES(x) 

        Parameters
        ----------
        x: integer
            ADC counts
        ll: float
            Expected value of the poisson law (avergae thermal photon number in 60ns)
        xt : float
            Cross talk probability
        g: float
            gain : Averaged number of ADC counts produced by one photo-electron
        sigma_pe: float
            st.d. of the ADC counts produced by one photo-electron
        sigma_el: float
            st.d. of the ADC counts produced by the electronic noise
        """

        STP = np.sqrt(2 * np.pi)

        true_Ext = 1/(1-xt) - ll /self.n_sample * self.peak_search_window_width
        x00 = -ll * g / self.n_sample * self.peak_search_window_width * true_Ext


        def single_gauss(x, xn, sigma_n):
            return np.exp(-1/2 *((x-xn)/sigma_n)**2) / (sigma_n * STP)

        def sigma_n(n, sigma_pe, sigma_el):
            return np.sqrt(n * sigma_pe**2 + sigma_el**2)
        def G_poisson(n,theta,l_lambda): ## not folowing notation from On the Generalized Poisson Distribution (2006) lambda = muXT!
            a1 = theta*(theta + n*l_lambda)**(n-1) / scsp.factorial(n)
            a2 = np.exp(-theta - n*l_lambda)
            return a1*a2

        def borel(k,m):
            return np.exp(-m*k)*(m*k)**(k-1) / scsp.factorial(k)



        def spe_spectrum_function_(x, ll, xt, g, x_el, sigma_pe, sigma_el):

            S = np.zeros_like(x)
            for n in range(0, 8):
                if n==0:
                    S+= single_gauss(x, x_el, sigma_el) * scst.poisson(ll).pmf(0)
                else:
                    xn =  x00 + n * g
                    _sigma_n = sigma_n(n, sigma_pe, sigma_el)
                    if True:
                        S += (borel(n,xt)*
                              (1-scst.poisson(ll).pmf(0))*
                              single_gauss(x, xn, _sigma_n)
                              )
                    if False:
                        S += (
                            G_poisson(n,ll,xt)*
                            single_gauss(x, xn, _sigma_n)
                            )

            return  S

        return spe_spectrum_function_(x, ll, xt, g, x_el, sigma_pe, sigma_el) * self.tot_evts* self.binwidth

    ##############################################

    ##############################################
    ##############################################
    ##############################################

    def spe_spectrum_function(self,x, ll, xt, g, sigma_pe, sigma_el,x_el):
        """
        single phot-electron spectrum fuction x -> MES(x) 

        Parameters
        ----------
        x: integer
            ADC counts
        ll: float
            Expected value of the poisson law (avergae thermal photon number in 60ns)
        xt : float
            Cross talk probability
        g: float
            gain : Averaged number of ADC counts produced by one photo-electron
        sigma_pe: float
            st.d. of the ADC counts produced by one photo-electron
        sigma_el: float
            st.d. of the ADC counts produced by the electronic noise
        """

        STP = np.sqrt(2 * np.pi)

        #true_Ext = 1/(1-xt) - ll /self.n_sample * self.peak_search_window_width
        x00 = -ll * g / self.n_sample * self.peak_search_window_width * 1/(1-xt)
        ## TNT
        ##x_el = sigma_el*2.17 ??

        def single_gauss(x, xn, sigma_n):
            return np.exp(-1/2 *((x-xn)/sigma_n)**2) / (sigma_n * STP)

        def sigma_n(n, sigma_pe, sigma_el):
            return np.sqrt(n * sigma_pe**2 + sigma_el**2)

        def G_poisson(n,theta,l_lambda): ##  folowing notation from On the Generalized Poisson Distribution (2006) lambda = muXT!
            a1 = theta*(theta + n*l_lambda)**(n-1) / scsp.factorial(n)
            a2 = np.exp(-theta - n*l_lambda)
            return a1*a2

        def borel(k,m):
            return np.exp(-m*k)*(m*k)**(k-1) / scsp.factorial(k)



        def spe_spectrum_function_(x, ll, xt, g, x_el, sigma_pe, sigma_el):
            poissonratio= (1-np.exp(-ll))/(1-np.exp(-ll/self.n_sample * self.peak_search_window_width))
            S = np.zeros_like(x)
            for n in range(0, 8):
                if n==0:
                    S+= single_gauss(x, x_el, sigma_el) * scst.poisson(ll).pmf(0)
                else:
                    xn =  x00 + n * g
                    _sigma_n = sigma_n(n, sigma_pe, sigma_el)
                    if False:
                        S += (borel(n,xt)*
                              (1-scst.poisson(ll).pmf(0))*
                              single_gauss(x, xn, _sigma_n)
                              )
                    if True:
                        S += (
                            G_poisson(n,ll/ self.n_sample * self.peak_search_window_width,xt)*
                            poissonratio*
                            single_gauss(x, xn, _sigma_n)
                            )

            return  S

        return spe_spectrum_function_(x, ll, xt, g, x_el, sigma_pe, sigma_el) * self.tot_evts* self.binwidth

    ##############################################
    def plot_failed_pixs(self):
        for pix in fitter.res['pixel_id'][fitter.res['calib_flag']==0] :
           fitter.plot_onepix(pix)

    ##############################################

    def aspe_fit(self,pix,min_ADC = -15):
        X = self.centers_adcsum[self.centers_adcsum > min_ADC]
        Y = self.Qsum_hist[pix][self.centers_adcsum > min_ADC]
        if Y[X>50].sum() >0:
            g0 = 20
        else:
            Y = Y[X<50]
            X = X[X<50]
            g0=0

        # X = X[Y>0]
        # Y = Y[Y>0]
        Yerr  = np.sqrt(Y)
        Yerr[Yerr==0] = 1
        # least_squares = LeastSquares(X, Y, Yerr, self.spe_spectrum_function)
        def likelihood(logL       ,
                       xt,
                       g       ,
                       sigma_pe,
                       sigma_el,
                       x_el
                       ):
            preds = self.spe_spectrum_function(X,logL,xt,g,sigma_pe,sigma_el,x_el)
            # l =  np.sum( [np.log(scipy.stats.poisson(preds[ii]).pmf(Y[ii])) for ii in range(len(Y)) ] )
            if ((preds>0).all()) & (np.isreal(preds).all()):
                logL =  np.sum( np.log(preds)*Y-preds-scipy.special.gammaln(Y+1) )
            else:
                logL=-10000
            return -2*logL

        m = Minuit(likelihood,
                   l        = .2,      # lambda
                   xt       = 0.11,      # X-talk
                   g        = g0,       # g
                   sigma_pe = 2 ,      # sigma_pe
                   sigma_el = 4.,       # sigma_el
                   x_el = 0,
                   # throw_nan=True
                   )
        #m.fixed['sigma_pe']=True
        #m.fixed['x_el']=False

        ## new version of Minuit :
        if True:
            m.limits['l']  = (0,2)
            m.limits['xt']  = (1e-5,0.4)
            m.limits['g']  = (0,40)
            m.limits['sigma_pe'] = (1,10)
            m.limits['sigma_el'] = (0,20)
        m.errordef=0.5   ### 0.5 for likelihood, 1 for LS
        m.migrad()
        m.hesse()






        # if m.accurate:
        #     result = dict(zip(self.param_names, m.np_values() ))
        #     result['P_chi2'] = (((Y - self.spe_spectrum_function(X, *m.np_values())) / Yerr )**2).sum() / (len(X) - len(self.param_names))
        # else:
        #     print('fit_failed')
        #     result = dict(zip(self.param_names, [-1,-1,-1,-1,-1]))
        #     result['P_chi2'] = 100


        result = dict()
        result['pixel_id']     = pix
        for ii, param in enumerate(self.param_names):
            result[param]        = m.values[ii]
            result[param+'_err'] = m.errors[m.parameters[ii]]
        ## reduced Chi2
        result['P_chi2'] = (((Y - self.spe_spectrum_function(X, *np.array(m.values) )) / Yerr )**2).sum() / (len(X) - len(self.param_names))



        true_Ext = 1/(1-result["X_talk"]) - result["dcr"]/self.n_sample * self.peak_search_window_width





        result['dcr_MHz']      = (m.values['l']*5) # MHz
        result['true_Xtalk']   = 1-1/true_Ext
        result['fit_accuracy'] = m.accurate
        result['dc_to_pe']     = result['gain']*1/(1-result["true_Xtalk"])
        result['x0']           = -result['dcr'] * result['dc_to_pe']/ self.n_sample * self.peak_search_window_width

        conf_limits = dict()
        conf_limits['dcr']      = [0,1.5]
        conf_limits['X_talk']   = [0.,0.3]
        conf_limits['gain']     = [16,35]
        conf_limits['sigma_pe'] = [1,5.5]
        conf_limits['sigma_el'] = [0,7]
        conf_limits['P_chi2']   = [0,15]
        conf_limits['true_Xtalk']   = [0.05,0.25]

        ## estimate the validity of fitted parameter
        ## (pixels with important electronic noise lead to degenerated shapes)
        calib_flag = int(m.accurate)
        for key in conf_limits.keys():
            if result[key]<conf_limits[key][0] or result[key]>conf_limits[key][1]:
                calib_flag = 0
        result['calib_flag']   = calib_flag

        self.minuit = m ## will not be filled in pool mp
        return result

    ##############################################

    def do_all_fit(self,n_proc=1):

        centers_adcsum,Qsum_hist,centers_adcmax,Qmax_hist = self.get_histograms()
        pool = mp.Pool(n_proc)
        pool_results = pool.map(self.aspe_fit,
                                self.pixels,
                                )
        results = dict()
        for key in pool_results[0].keys():
            results[key] = []

        for result in pool_results:
            for key in result.keys():
                    results[key].append(result[key])

        self.res            = pd.DataFrame(results)

        self.results = results
        pool.close()
        self.n_failed = (self.res['calib_flag']<1).sum()
        print('Calibration Done!')
        print(f'Calibration failed for {self.n_failed}  pixels ({self.n_failed/self.n_pixels*100:.3} %)')

        return

    ##############################################


    def plot_cam_dist(self,save_plots=False,nan_failed=False):

        for key in self.results.keys():
            try:
                f,ax = plt.subplots()
                disp = CameraDisplay(self.geom,ax=ax)
                disp.add_colorbar()
                image = self.res[key].copy()
                if nan_failed:
                    image[self.res['calib_flag']==0] = np.nan
                disp.image = image
                ax.set_title(key+f" Tel {self.tel}")
                #disp.set_limits_percent(95)
                disp.set_limits_minmax(image.min()*0.95,image.max()*1.05)
                if save_plots:
                    f.savefig(self.plot_dir+f'cam_{key}_tel{self.data_path[-2]}.png')
            except Exception:
                print("Failed")
                return
    ##############################################

    def plot_onepix(self,pix): ## todo
        result = dict(self.res.iloc[pix])

        # fitlabel = 'fit : \n gain : {:.2} \n Xt : {:.2}'.format(result['gain'],
        #                                                         result['B_param'])
        fitlabel = 'fit : \n'
        for key in self.param_names:
            fitlabel=fitlabel+'{}    ::   {:.3} $\\pm {:.2}$ \n'.format(key,result[key],result[key+'_err'])
        fitlabel = fitlabel+' Calib Flag :: {}'.format(result["calib_flag"])
        print(self.res.iloc[pix])
        f,ax = plt.subplots(figsize=(12,5))
        # ax.plot(self.centers_adcsum,self.Qsum_hist[pix],label='$\sum$ ADC - baseline')
        ax.fill_between(self.centers_adcsum,
                        self.Qsum_hist[pix]+(self.Qsum_hist[pix])**0.5,
                        self.Qsum_hist[pix]-(self.Qsum_hist[pix])**0.5,
                        alpha = .5,
                        color = 'green',
                        label=r'$\sum$ ADC - baseline')

        ax.plot(self.centers_adcsum,self.spe_spectrum_function(self.centers_adcsum,
                                                            *[result[k] for k in self.param_names]),
                                                            '--',
                                                            color='black',
                                                            label= fitlabel)

        def single_gauss(x,A,x0,sigma_el):
            return np.exp(-1/2 *((x-x0)/sigma_el)**2) / (sigma_el * np.sqrt(2 * np.pi)) *A


        ax.plot(self.centers_adcsum,
                single_gauss(self.centers_adcsum,
                             self.tot_evts* np.exp(-result['dcr']),
                             result['x_el'],
                             result['sigma_el']),
                ':',
                label='pedestal')

        ax.set_yscale('log')
        ax.set_ylim(1e-1,self.Qsum_hist[pix].max()*2)
        ax.grid()
        ax.legend()
        ax.set_xlabel(r'$\Sigma$ ADC')
        ax.set_title(f'SPE spectrum -- Tel {self.tel} -- pix {pix}')

        return f,ax

    ##############################################

    def plot_param_dist(self,save_plot=False,show_uncal=True):
        ## plot hist of param distribution
        mask = self.res['calib_flag']==1

        kwargs = dict(histtype='stepfilled', alpha=0.3, ec="k")
        # kwargs = dict(histtype='step',alpha=.8)
        for key in self.res.keys()[:-1]:
            f,ax = plt.subplots()
            bins = np.linspace(self.res[key].min(),
                           self.res[key].max(),
                           100)

            ax.hist(self.res[key][mask],
                bins=bins,
                label = f'tel2 -- median = {np.median(self.res[key]):.3}',
                **kwargs)
            if show_uncal:
                ax.hist(self.res[key][~mask],
                    bins=bins,
                    color='red',
                    **kwargs)

            ax.set_title(key)
            ax.legend(loc='upper left')
            ax.grid()
            if key not in[ 'B_param', 'gain','P_chi2']:
                ax.set_xlabel('ADC count')
            if key =='gain':
                ax.set_xlabel('ADC / p.e.')
            if save_plot:
                f.savefig(self.plot_dir+f'hist_{key}.png')

    ##############################################

    def plot_all_ses(self):
        ### plot ses for all pix
        f, ax = plt.subplots(figsize=(12,6))
        for pix in range(self.Qsum_hist.shape[0]):
            ax.plot(self.centers_adcsum,
                    self.Qsum_hist[pix],
                    color='black',
                    alpha=.1)
        ax.grid()
        ax.set_xlabel(r'$\Sigma$ ADC')
        ax.set_yscale('log')

    ##############################################


    def save_h5(self,res_save_name=None):
        if res_save_name is None:
            res_save_name = self.save_str

        df_param = pd.DataFrame.from_dict(self.results)



        df_hists = pd.DataFrame(self.Qsum_hist[:,:])
        df_hists.columns = self.centers_adcsum

        df_param.to_hdf(os.path.join(self.save_dir,
                                     res_save_name+'fitted_parameters.h5'),
                        "df_param")
        df_hists.to_hdf(os.path.join(self.save_dir,
                                     res_save_name+'histograms.h5'),
                        'df_hists')

        return

    # def read_h5(self,h5file_basename):
    #     h5filepath = os.join(self.save_dir,h5file_basename)
    #     df_hists = pd.read_hdf(h5filepath+'fitted_parameters.h5')
    #     df_param = pd.read_hdf(h5filepath+'histograms.h5')

    #     self.results = df_param.to_dict(orient='list')
    #     # self.centers_adcsum = np.array(df_hists.indexes)
    #     self.Qsum_hist = np.array(df_hists)


#########################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--plot_dir',   type=str , default = './SPE_plots_mc/' )
    parser.add_argument('--data_path',  type=str , default = '/net/' )
    parser.add_argument('--day',        type=int , default = 2)
    parser.add_argument('--month',      type=int , default = 3)
    parser.add_argument('--year',       type=int , default = 2023)
    parser.add_argument('--first_file', type=int , default = 136)
    parser.add_argument('--tel',        type=int , default = 1)
    parser.add_argument('--n_file',     type=int , default = 2)
    parser.add_argument('--save_dir',   type=str , default = './calib_h5/')
    parser.add_argument('--save_plot',  type=bool, default = False)
    parser.add_argument('--n_proc',     type=int , default = 20)

    args = parser.parse_args()



    fitter= mes_fitter( day          = args.day,
                        month        = args.month,
                        year         = args.year,
                        data_path    = args.data_path,
                        tel          = args.tel,
                        first_file_n = args.first_file,
                        n_files      = args.n_file,
                        max_evt      = 100000,
                        plot_dir     = args.plot_dir)

    r = fitter.do_all_fit(n_proc=args.n_proc)

    if args.save_dir is not None:
        fitter.save_h5()


