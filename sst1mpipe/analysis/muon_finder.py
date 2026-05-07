import argparse
import datetime
import glob
import multiprocessing as mp
import os
from os import path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pkg_resources
from ctapipe.image import tailcuts_clean
from ctapipe.image.cleaning import number_of_islands
from ctapipe.image.muon import ring_completeness
from ctapipe.instrument import SubarrayDescription
from ctapipe.io import EventSource
from ctapipe.visualization import CameraDisplay
from scipy.ndimage import convolve1d
from scipy.optimize import minimize

from sst1mpipe.calib.calib import get_default_window
from sst1mpipe.io import load_config
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.utils.monitoring_pedestals import sliding_pedestals
from sst1mpipe.utils.NSB_tools import VAR_to_Idrop

DEFAULT_CONFIG_FILE = pkg_resources.resource_filename(
                            'sst1mpipe',
                            path.join(
                                'data',
                                'sst1mpipe_data_config.json'
                            )
                        )
DEFAULT_CONFIG = load_config(DEFAULT_CONFIG_FILE, ismc=False)
def swap_r0_modules_59_88(event, tel=None):

    # module 59
    mask59 = np.zeros(1296, dtype=bool)
    mask59[1029] = True
    mask59[1098:1102+1] = True
    mask59[1133:1134+1] = True
    mask59[1064:1067+1] = True
    waveform_59 = event.sst1m.r0.tel[tel].adc_samples[mask59, :]
    bls59 = event.sst1m.r0.tel[tel].digicam_baseline[mask59]

    # module 88
    mask88 = np.zeros(1296, dtype=bool)
    mask88[1103] = True
    mask88[1165:1169+1] = True
    mask88[1194:1195+1] = True
    mask88[1135:1138+1] = True
    waveform_88 = event.sst1m.r0.tel[tel].adc_samples[mask88, :]
    bls88 = event.sst1m.r0.tel[tel].digicam_baseline[mask88]

    event.sst1m.r0.tel[tel].adc_samples[mask59] = waveform_88
    event.sst1m.r0.tel[tel].adc_samples[mask88] = waveform_59

    event.sst1m.r0.tel[tel].digicam_baseline[mask59] = bls88
    event.sst1m.r0.tel[tel].digicam_baseline[mask88] = bls59

    return event

def get_res_dict(ismc=False):
    result_dict  = dict({ 'toa'       : [],
                          'x'         : [],
                          'y'         : [],
                          'n_fmask'   : [],
                          'radius'    : [],
                          'mu_pixs'   : [],
                          'Qs'        : [],
                          'phi_pixs'  : [],
                          '5pe_ratio' : [],
                          'Q_mu'      : [],
                          'Q_nmu'     : [],
                          'ring_comp' : [],
                          'WM_pixs'   : [],
                          'WM_Qs'     : [],
                          'Q_raw'     : [],
                          'tom_std'   : [],
                          'event_id'  : [],
                          })
    if ismc:
        result_dict['E_mu']     = []
    else:
        result_dict['mbs']   = []
        result_dict['bsstd'] = []
    return result_dict

class muon_finder:
    def __init__(self,
                 filename    = '/net/cs2/data/raw/2024/01/29/SST1M2/SST1M2_20240129_0050.fits.fz',
                 plot_dir     = './muons_plots/',
                 tel          = 1,
                 max_evt      = 200,
                 gain         = 1.,
                 plot         = False,
                 ismc         = False):




        self.tailcuts_1 = 5
        self.tailcuts_2 = 4
        self.min_pix = 13
        self.max_pix = 50

        self.filename=filename
        if ismc:
            self.date_str =  filename.split('_')[-1].split(".")[0]
            self.file_n   =  int(filename.split('_')[-1].split(".")[0])
            self.day             = 0
            self.month           = 0
            self.year            = 0
        else:
            self.date_str =  filename.split('_')[-2]
            self.file_n = int(filename.split('_')[-1].split(".")[0])
            self.day             = int(self.date_str[6:])
            self.month           = int(self.date_str[4:6])
            self.year            = int(self.date_str[:-4])
            self.window_t,_      = get_default_window(tel)

        self.plot_dir        = plot_dir
        self.max_evt         = max_evt
        self.gain            = gain
        self.tel             = tel

        self.plot=plot
        self.ismc = ismc
        if self.plot and (not os.path.exists(self.plot_dir)):
            os.makedirs(self.plot_dir)

        if self.ismc :
            self.w_start = 16
            self.w_end   = 25
        else :
            self.w_start = 6
            self.w_end   = 15

        self.mu_data = get_res_dict(self.ismc)




        subarray_file = pkg_resources.resource_filename(
                            'sst1mpipe',
                            path.join(
                                'data',
                                'sst1m_array.h5'
                            )
                        )

        subarray  = SubarrayDescription.from_hdf(subarray_file, focal_length_choice="EQUIVALENT")
        self.geom = subarray.camera_types[0].geometry
        self.n_pixels = self.geom.n_pixels
        self.pixels   = np.arange(self.n_pixels)

        self.file_list = []

        self.mbs   = []
        self.bsstd = []




    def get_barytom(self,wf):

        Q = np.array([ wf[ii:ii+3].sum() for ii in range(len(wf)-3)])
        im = Q.argmax()
        return (wf[im:im+3] * np.arange(im,im+3)/wf[im:im+3].max() ).mean()

    def get_muons(self):
        def get_circle_dists(C,pixs):
            Cx,Cy,Cr = C
            dists = np.array([])
            for pix in pixs:
                Cc = np.array([Cx,Cy])
                P  = np.array([self.geom.pix_x[pix].value,
                               self.geom.pix_y[pix].value])
                d = np.abs(np.linalg.norm(Cc-P)-Cr)
                dists = np.append(dists,d)
            return dists.sum()

        def get_Weighted_circle_dists(C,pixs,Qs):
            Cx,Cy,Cr = C
            wdists = np.array([])
            for ii,pix in enumerate(pixs):
                Cc = np.array([Cx,Cy])
                P  = np.array([self.geom.pix_x[pix].value,
                               self.geom.pix_y[pix].value])
                d = np.abs(np.linalg.norm(Cc-P)-Cr)
                wdists = np.append(wdists,d*Qs.clip(0,30)[ii])
            return wdists.sum()

        if self.ismc:
            data_stream = EventSource(self.filename,max_events=self.max_evt)
            print(f"file {self.filename} opened")
        else :
            data_stream = SST1MEventSource(
                filelist    = [self.filename],
                max_events  = self.max_evt,
                disable_bar = True
                )

            pedestal_info = sliding_pedestals(input_file = self.filename,config=DEFAULT_CONFIG)
            pedestal_info.load_firsts_pedestals()
            if pedestal_info.get_n_events() == 0:
                print("No pedestal events found in firsts events. Skipping run")
                return
            else:
                print(f"{pedestal_info.get_n_events()} pedestals events loaded in buffer")

        #####################################
        ### Loop in all events to find muons:
        #####################################

        for ii,event in enumerate(data_stream):
            if self.ismc:
                if self.tel not in event.trigger.tels_with_trigger:
                    continue
                r0data = event.r0.tel[self.tel]
                r1data = event.r1.tel[self.tel]
                Q_sum_ADC    = (r0data.waveform.T[self.w_start:self.w_end]).sum(axis=0)
                Q_sum_window = (r1data.waveform.T[self.w_start:self.w_end]).sum(axis=0)
                E_mu = event.simulation.shower.energy
            else:


                tel = event.sst1m.r0.tels_with_data[0]
                if ii==0:
                    T0 = event.sst1m.r0.tel[tel].local_camera_clock/1e9
                    start_date = datetime.datetime.fromtimestamp(T0)
                    # datestr = "{}/{}/{} at {}h{}".format(start_date.day,
                    #                                      start_date.month,
                    #                                      start_date.year,
                    #                                      start_date.hour,
                    #                                      start_date.minute)
                    night_datestr = f"{start_date.day}-{start_date.day+1}/{start_date.month}/{start_date.year}"


                    print("night : "+night_datestr)
                if (tel==22) and (start_date<datetime.datetime(2024,7,18)):
                    event = swap_r0_modules_59_88(event, tel=tel)
                    #pass
                r0data = event.sst1m.r0.tel[tel]

                if r0data._camera_event_type.value==8:
                    pedestal_info.add_ped_evt(event, store_image=False)
                    pedestal_info.fill_mon_container(event)

                    mbs   = r0data.adc_samples.mean(axis=1)
                    bsstd = r0data.adc_samples.std(axis=1)
                    self.mbs.append(mbs[mbs>0].mean())
                    self.bsstd.append(bsstd[bsstd>2].mean())
                    if len(self.mbs)>100:
                        self.mbs   = self.mbs[-100:]
                        self.bsstd = self.bsstd[-100:]
                    continue
                else:
                    ## intergrate signal in a fixed window :
                    VI = VAR_to_Idrop(pedestal_info.get_charge_std().mean()**2, 20+self.tel)
                    Q_sum_ADC    = (r0data.adc_samples.T[self.w_start:self.w_end] - r0data.digicam_baseline).sum(axis=0)
                    Q_sum_window = Q_sum_ADC /self.gain /VI /self.window_t


            mask_tailcuts = tailcuts_clean(self.geom,
                                           Q_sum_window,
                                           self.tailcuts_1,
                                           self.tailcuts_2,
                                           min_number_picture_neighbors=2)


            n_island,mask_island = number_of_islands(self.geom, mask_tailcuts)

            #drop empty
            if n_island == 0:
                continue

            final_mask = mask_tailcuts



            cleaned_image = Q_sum_window.copy()
            cleaned_image[~final_mask] = 0



            if final_mask.sum() > self.min_pix and final_mask.sum() < self.max_pix:


                x0 = [self.geom.pix_x[final_mask].value.mean(),
                      self.geom.pix_y[final_mask].value.mean(),
                      .1]

                res = minimize(get_circle_dists, x0,args=(self.pixels[final_mask]), method='BFGS', tol=1e-2)
                #print(res.x)


                # rc = ring_completeness( self.geom.pix_x.value,
                #                         self.geom.pix_x.value,
                #                         weights  = Q_sum_window,
                #                         radius   = res.x[2],
                #                         center_x = res.x[0],
                #                         center_y = res.x[1],
                #                         threshold=4,
                #                         bins=30
                #                         )
                # print(ii,rc)


                # if res.success:
                # if res.fun<0.2:
                # if res.x[2]>0.1:
                # if rc>0.4:
                ## Simple preselection
                ismu_test = (res.x[2]>0.05) & (res.x[2]<.15)
                ismu_test&= Q_sum_window.max() < 60.

                if ismu_test:
                    Cx,Cy,Cr = res.x

                    d_center = ( (self.geom.pix_x.value-Cx)**2 +
                                 (self.geom.pix_y.value-Cy)**2)**.5

                    WM_dr = 0.04
                    wide_mask = (d_center > Cr-WM_dr) & (d_center < Cr+WM_dr)

                    ## refit the circle
                    res = minimize(get_Weighted_circle_dists,
                                   res.x,
                                   args=(self.pixels[wide_mask],Q_sum_window[wide_mask]),
                                   method='BFGS',
                                   tol=1e-2)
                    Cx,Cy,Cr = res.x



                    d_center = ( (self.geom.pix_x.value-Cx)**2 +
                                 (self.geom.pix_y.value-Cy)**2)**.5


                    dr = 0.015
                    mu_mask = (d_center > Cr-dr) & (d_center < Cr+dr)

                    ## if there is signal outside the ring : skip the event
                    if Q_sum_window[~mu_mask & final_mask].sum() > 60 :
                        continue
                    rc = ring_completeness(self.geom.pix_x.value[mu_mask],
                                           self.geom.pix_y.value[mu_mask],
                                           weights  = Q_sum_window[mu_mask],
                                           radius   = res.x[2],
                                           center_x = res.x[0],
                                           center_y = res.x[1],
                                           threshold=7,
                                           bins=12
                                           )
                    ## if there is less than half a muon : skip the event
                    if rc < 0.5:
                        continue

                    ## looking at the time dispersion.. Is it useful? If not, I'll remove this
                    if False:
                        wfs  = r0data.adc_samples[mu_mask * Q_sum_window > 8, 4:16]
                        peak_search_window_width = 3
                        sums = convolve1d(
                            wfs, np.ones(peak_search_window_width), axis=1, mode="nearest"
                        )
                        peak_index = np.argmax(sums[:, 2:-2], axis=1) + 2
                        bin_time = 4
                        rawtoms = np.array([
                            np.sum(
                                wfs[ii, m - 1 : m + 2]
                                / wfs[ii, m - 1 : m + 2].sum()
                                * np.arange((m - 1) * bin_time, (m + 2) * bin_time, bin_time)
                            )
                            if wfs[ii, m - 1 : m + 2].sum() > 0
                            else m*bin_time
                            for ii, m in enumerate(peak_index)
                            ])
                        if rawtoms.shape[0]>2:
                            tom_std = rawtoms.std()
                        else:
                            tom_std = 0
                    tom_std = 0

                    phi_pixs = [np.arctan2(self.geom.pix_y[p].to('m').value-Cy,
                                           self.geom.pix_x[p].to('m').value-Cx) for p in self.pixels[mu_mask]]


                    self.mu_data['x'].append(Cx)
                    self.mu_data['y'].append(Cy)
                    self.mu_data['n_fmask'].append(final_mask.sum())
                    self.mu_data['radius'].append(Cr)
                    self.mu_data['mu_pixs'].append(self.pixels[mu_mask])
                    self.mu_data['Qs'].append(Q_sum_window[mu_mask])
                    self.mu_data['phi_pixs'].append(np.array(phi_pixs))
                    self.mu_data['5pe_ratio'].append( (Q_sum_window[mu_mask] > 5).sum() / mu_mask.sum() )
                    self.mu_data['Q_mu'].append(Q_sum_window[mu_mask].sum())
                    self.mu_data['Q_raw'].append(Q_sum_ADC[mu_mask].sum())
                    self.mu_data['Q_nmu'].append(Q_sum_window[~mu_mask & final_mask].sum())
                    self.mu_data['WM_pixs'].append(self.pixels[wide_mask])
                    self.mu_data['WM_Qs'].append(Q_sum_window[wide_mask])
                    self.mu_data['tom_std'].append(tom_std)
                    self.mu_data['ring_comp'].append(rc)
                    if self.ismc :
                        self.mu_data['E_mu'].append(E_mu.to('GeV').value)
                        self.mu_data['event_id'].append(ii)
                        self.mu_data['toa'].append(ii)
                    else:
                        self.mu_data['mbs'].append(np.median(self.mbs))
                        self.mu_data['bsstd'].append(np.median(self.bsstd))
                        self.mu_data['event_id'].append(event.sst1m.r0.event_id)
                        self.mu_data['toa'].append(r0data.local_camera_clock/1e9)

                    #self.mbs   = []
                    #self.bsstd = []

                    for key in self.mu_data.keys():
                        try:
                            print(key+str(self.mu_data[key][-1]))
                        except (IndexError, TypeError):
                            pass
                    if (self.plot) and (np.random.randint(100)<110) and (np.median(self.bsstd)>11) and (rc>0.6):
                        f,axs = self.plot_event(Q_sum_window,final_mask,mu_mask)
                        f.savefig(self.plot_dir+f'd{self.date_str}_f{self.file_n}_ev_{ii}_ismuontest_rc{rc:.4}.png' )
                        plt.close(f)

                    if  len(self.mu_data['event_id'])%500==0 and False:
                        df = pd.DataFrame.from_dict(self.mu_data)
                        df.to_hdf(os.path.join('./',
                                               f'mudata_tel{self.tel}_t.h5'),
                                               "df")



        return

    def plot_event(self,Q_sum,mask,mu_mask):
        f,axs = plt.subplots(1,3,figsize=(15,4))
        # disp = CameraDisplay(digicam_geometry,ax=ax)
        disp = CameraDisplay(self.geom,ax=axs[0])
        disp.image = Q_sum
        disp.add_colorbar(ax=axs[0])

        image = np.zeros(self.n_pixels)
        image[mask]=Q_sum[mask]
        disp = CameraDisplay(self.geom,ax=axs[1])
        disp.image = image
        axs[1].set_title(f"Qtot = {Q_sum.sum()} \n Qmax = {Q_sum.max()} (b std ={np.median(self.bsstd):.3})")
        disp.add_colorbar(ax=axs[1])


        image = np.zeros(self.n_pixels)
        image[mu_mask]=Q_sum[mu_mask]
        disp = CameraDisplay(self.geom,ax=axs[2])

        disp.image =image
        disp.highlight_pixels(self.pixels[mu_mask],'red',.5)
        disp.add_colorbar(ax=axs[2])

        ttt = ""
        for key in ['radius','5pe_ratio','Q_mu','Q_nmu']:
            ttt = ttt+(f"{key} : {self.mu_data[key][-1]:.4} \n" )

        tetas = np.linspace(0,2*np.pi,30)
        Cx,Cy,Cr = [self.mu_data[key][-1] for key in ['x','y','radius']]
        # axs[0].plot(Cx+np.cos(tetas)*Cr,Cy+np.sin(tetas)*Cr,color='grey')
        axs[1].plot(Cx+np.cos(tetas)*Cr,Cy+np.sin(tetas)*Cr,color='grey',label=ttt)

        # axs[2].set_title(ttt)
        axs[1].legend()
        return f,axs








if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    if True:
        parser.add_argument('--data_path',    type=str, default = "")
        parser.add_argument('--plot_dir',     type=str, default = None )
        parser.add_argument('--table_dir',    type=str, default = '/mnt/nfs28_auger3/tavernier/sst1m/muon_ana/muons_table20241126' )
        parser.add_argument('--tel',          type=int, default = 2)
        parser.add_argument('--max_evt',      type=int, default = 20000000)
        parser.add_argument('--ismc',         type=bool, default = False)


    args = parser.parse_args()


    if args.tel==2:
        calib_df = pd.read_hdf('/mnt/nfs19/tavernier/sst1m/calib/averaged_calib_param_v2_tel2.h5')
    else:
        calib_df = pd.read_hdf('/mnt/nfs19/tavernier/sst1m/calib/averaged_calib_param_v2_tel1.h5')
    gains = np.array(calib_df['dc_to_pe'])
    gains[~np.isfinite(gains)]=gains[np.isfinite(gains)].mean()

    def proces_one_file(filename):
        finder = muon_finder( filename     = filename,
                              plot_dir     = args.plot_dir,
                              tel          = args.tel,
                              max_evt      = args.max_evt,
                              ismc         = args.ismc,
                              gain         = gains)

        try:
            finder.get_muons()
        except Exception:
            print(f"file {filename} FAILED")


        return finder.mu_data

    n_proc = 8
    print("ready")

    results = get_res_dict(args.ismc)

    file_list=glob.glob(args.data_path+'/*.?z')

    pool = mp.Pool(n_proc)
    pool_results = pool.map(proces_one_file,
                            file_list[:1500])
    pool.close()



    for res in pool_results:
        for key in results.keys():
            results[key] = results[key]+res[key]



    # hdf5_file_path = os.path.join(args.table_dir,"tmp_"+args.out_file)
    # df.to_hdf(hdf5_file_path, key='df', mode='w')
    # df = pd.DataFrame(dict(results)).sort_values(by=['time'])
    df = pd.DataFrame(dict(results)).sort_values(by=['toa'])
    if args.ismc:
        df.to_hdf(os.path.join(f'{args.table_dir}',
                               'MC_{}_{}_mudataV30_tel{}.h5'.format(args.data_path.split("/")[-3],
                                                                  args.data_path.split("/")[-2],
                                                                  args.tel)),
                               "df")
    else:
        df.to_hdf(os.path.join(f'{args.table_dir}',
                               '{}{}{}_mudataV30_tel{}.h5'.format(args.data_path.split("/")[7],
                                                                  args.data_path.split("/")[8],
                                                                  args.data_path.split("/")[9],
                                                                  args.tel)),
                               "df")
