import logging
from collections import deque
from statistics import mean

import astropy.units as u
import numpy as np
from ctapipe.calib import CameraCalibrator
from ctapipe.image import ImageProcessor

from sst1mpipe.calib import (
    Calibrator_R0_R1,
    get_window_corr_factors,
    saturated_charge_correction,
    window_transmittance_correction,
)
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.utils import get_subarray, get_swap_flag

MON_EVT_TYPE = 8
class sliding_pedestals:
    def __init__(self, input_file=None, max_array_size = 100, max_images_array = 1000, config=None):

        self.timestamps = deque([])
        self.ped_mean_array   = deque([])
        self.ped_std_array    = deque([])
        self.ped_img_array    = deque([])
        self.max_array_size   = max_array_size
        self.max_images_array = max_images_array
        self.ped_img_sum      = None
        self.ped_img_sum2     = None
        self.processed_pedestals = 0
        self.config = config
        self.input_file = input_file
        self.pedestals_in_file = True

        self.log_pedestal_settings()
        self.load_firsts_pedestals()

        if self.get_n_events() == 0:
            logging.warning("No pedestal events found in firsts events. Cleaned shower/NSB events used instead.")
            self.load_firsts_fake_pedestals()
            logging.info(f"{self.get_n_events()} fake pedestals events loaded in buffer")
            self.pedestals_in_file = False
        else:
            logging.info(f"{self.get_n_events()} pedestals events loaded in buffer")
            self.pedestals_in_file = True

    def add_ped_evt(self, evt, cleaning_mask=None, store_image=True):
        tel = evt.sst1m.r0.tels_with_data[0]
        pedestal = evt.sst1m.r0.tel[tel].adc_samples
        if store_image:
            image = evt.dl1.tel[tel].image
        else:
            image = np.zeros(pedestal.shape[0], dtype=np.float64)
        self.processed_pedestals = self.processed_pedestals + 1
        if cleaning_mask is not None:
            pedestal[cleaning_mask] = -100 * np.ones(pedestal.shape[1])
            image[cleaning_mask] = -100
        timestamp = evt.sst1m.r0.tel[tel].local_camera_clock/1e9

        self.ped_img_array.append(image)
        if self.ped_img_sum is None:
            self.ped_img_sum = image
            self.ped_img_sum2 = image**2
        else:
            self.ped_img_sum += image
            self.ped_img_sum2 += image**2

        if len(self.ped_img_array) > self.max_images_array:
            self.ped_img_sum -= self.ped_img_array[0]
            self.ped_img_sum2 -= self.ped_img_array[0]**2
            self.ped_img_array.popleft()

        self.ped_mean_array.append(pedestal.mean(axis=1))
        self.ped_std_array.append(pedestal.std(axis=1))
        self.timestamps.append(timestamp)
        if len(self.timestamps) > self.max_array_size:
            self.ped_mean_array.popleft()
            self.ped_std_array.popleft()
            self.timestamps.popleft()

    def get_n_events(self):
        return len(self.timestamps)

    def get_mean_ts(self):
        return mean(self.timestamps) *u.s

    def get_min_ts(self):
        return self.timestamps[0] *u.s

    def get_max_ts(self):
        return self.timestamps[-1] *u.s

    def get_charge_mean(self):
        pedarray = np.array(self.ped_mean_array)
        masked = np.ma.masked_values(pedarray, -100)
        return masked.mean(axis=0).data

    def get_charge_median(self):
        pedarray = np.array(self.ped_mean_array)
        masked = np.ma.masked_values(pedarray, -100)
        return np.median(masked.data,axis=0)

    def get_charge_std(self):
        pedarray = np.array(self.ped_std_array)
        masked = np.ma.masked_values(pedarray, 0)
        return masked.data.mean(axis=0)

    def get_img_charge_mean(self):
        return self.ped_img_sum/len(self.ped_img_array)

    def get_img_charge_std(self):
        return np.sqrt(self.ped_img_sum2/len(self.ped_img_array) - (self.ped_img_sum/len(self.ped_img_array))**2)

    def fill_mon_container(self, evt):
        tel = evt.sst1m.r0.tels_with_data[0]
        mon_container=evt.mon.tel[tel].pedestal
        mon_container.n_events        = self.get_n_events()
        mon_container.sample_time     = self.get_mean_ts()
        mon_container.sample_time_min = self.get_min_ts()
        mon_container.sample_time_max = self.get_max_ts()
        mon_container.charge_mean     = self.get_charge_mean()
        mon_container.charge_median   = self.get_charge_median()
        mon_container.charge_std      = self.get_charge_std()
        return

    def load_firsts_pedestals(self,max_n_ped=100,max_n_img=1000,max_evt=100000):
        """
        Reads first max_n_ped pedestal events in the buffer.
        """

        data_stream = SST1MEventSource([self.input_file],
                                       max_events=max_evt)

        for ii,event in enumerate(data_stream):

            if ii == 0:
                tel = event.sst1m.r0.tels_with_data[0]

            r0data = event.sst1m.r0.tel[tel]
            if r0data._camera_event_type.value == MON_EVT_TYPE:
                self.add_ped_evt(event, store_image=False)

            if len(self.timestamps) >= max_n_ped:
                break

        data_stream._subarray = get_subarray()
        r1_dl1_calibrator = CameraCalibrator(subarray=data_stream.subarray, config=self.config)

        # to treat images we need to estimate gain drop from pedestals, so redo pedestals once more...
        keep_size = self.max_array_size
        self.max_array_size = self.max_images_array
        for jj,event in enumerate(data_stream):

            if jj == 0:
                calibrator_r0_r1 = Calibrator_R0_R1(config=self.config, telescope=tel)
                window_corr_factors, _ = get_window_corr_factors(
                            telescope=tel, config=self.config
                            )
                swap_modules = get_swap_flag(event)

            r0data = event.sst1m.r0.tel[tel]
            if r0data._camera_event_type.value == MON_EVT_TYPE:

                # here we apply gain drop correction
                event = calibrator_r0_r1.calibrate(event, pedestal_info=self)
                event.r1.tel[tel].selected_gain_channel = np.zeros(data_stream.subarray.tel[tel].camera.readout.n_pixels,dtype='int8')

                r1_dl1_calibrator(event)

                # Integration correction of saturated pixels
                event = saturated_charge_correction(event)

                event = window_transmittance_correction(
                    event,
                    window_corr_factors=window_corr_factors,
                    telescope=tel,
                    swap_flag=swap_modules
                    )

                if jj <= ii:
                    self.ped_mean_array.popleft()
                    self.ped_std_array.popleft()
                    self.timestamps.popleft()
                    self.ped_img_array.popleft()

                self.add_ped_evt(event)

                if len(self.ped_img_array) >= max_n_img:
                    break

        self.max_array_size = keep_size
        for i in range(len(self.timestamps)-keep_size):
            self.ped_mean_array.pop()
            self.ped_std_array.pop()
            self.timestamps.pop()

    def load_firsts_fake_pedestals(self, max_evt=10):
        """
        Reads first max_evt fake pedestal events in the buffer.
        Fake pedestal is a shower event with Cherenkov pixels masked out.
        """

        source = SST1MEventSource([self.input_file],
                                       max_events=max_evt)
        source._subarray = get_subarray()
        r1_dl1_calibrator = CameraCalibrator(subarray=source.subarray, config=self.config)

        # Here (for the first few events) we use just the simple ImageProcessor, nothing fancy
        self.config["ImageProcessor"]["image_cleaner_type"] = "TailcutsImageCleaner"
        image_processor   = ImageProcessor(subarray=source.subarray, config=self.config)

        for ii,event in enumerate(source):

            if ii == 0:
                tel = event.sst1m.r0.tels_with_data[0]
                calibrator_r0_r1 = Calibrator_R0_R1(config=self.config, telescope=tel)

            event = calibrator_r0_r1.calibrate(event)
            event.r1.tel[tel].selected_gain_channel = np.zeros(source.subarray.tel[tel].camera.readout.n_pixels,dtype='int8')

            r1_dl1_calibrator(event)
            image_processor(event)

            clenaning_mask = event.dl1.tel[tel].image_mask
            # Arbitrary cut, just to prevent too big showers from being used 
            if sum(clenaning_mask) < 20:
                self.add_ped_evt(event, cleaning_mask=clenaning_mask, store_image=False)

        # to treat images we need to estimate gain drop from pedestals, so redo pedestals once more...
        keep_size = self.max_array_size
        self.max_array_size = self.max_images_array
        for jj,event in enumerate(source):

            if jj == 0:
                window_corr_factors, _ = get_window_corr_factors(
                            telescope=tel, config=self.config
                            )
                swap_modules = get_swap_flag(event)

            # here we apply gain drop correction
            event = calibrator_r0_r1.calibrate(event, pedestal_info=self)
            event.r1.tel[tel].selected_gain_channel = np.zeros(source.subarray.tel[tel].camera.readout.n_pixels,dtype='int8')

            r1_dl1_calibrator(event)
            image_processor(event)

            clenaning_mask = event.dl1.tel[tel].image_mask
            # Arbitrary cut, just to prevent too big showers from being used 
            if sum(clenaning_mask) < 20:
              # Integration correction of saturated pixels - done only here because the fake pedestals must match in both loops
              event = saturated_charge_correction(event)

              event = window_transmittance_correction(
                  event,
                  window_corr_factors=window_corr_factors,
                  telescope=tel,
                  swap_flag=swap_modules
                  )

              if jj <= ii:
                  self.ped_mean_array.popleft()
                  self.ped_std_array.popleft()
                  self.timestamps.popleft()
                  self.ped_img_array.popleft()

              self.add_ped_evt(event)

              if len(self.ped_img_array) >= self.max_images_array:
                  break

        self.max_array_size = keep_size
        for i in range(len(self.timestamps)-keep_size):
            self.ped_mean_array.pop()
            self.ped_std_array.pop()
            self.timestamps.pop()


    def log_pedestal_settings(self):

        if self.config['NsbCalibrator']['apply_pixelwise_Vdrop_correction']:
            logging.info("Voltage drop correction is applied pixelwise")

        if self.config['NsbCalibrator']['apply_global_Vdrop_correction']:
            logging.info("Voltage drop correction is applied globaly")

        if self.config['NsbCalibrator']['apply_global_Vdrop_correction'] == self.config['NsbCalibrator']['apply_pixelwise_Vdrop_correction']:
            if self.config['NsbCalibrator']['apply_global_Vdrop_correction']:
                logging.error("Voltage drop correction is applied 2 times!!! this is WRONG!")
            else:
                logging.warning("NO Voltage drop correction is applied")
