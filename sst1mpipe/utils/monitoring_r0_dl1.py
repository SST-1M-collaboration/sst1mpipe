import logging

import numpy as np

from sst1mpipe.utils import get_target


class Monitoring_R0_DL1:

    """
    Monitoring data for R0 -> DL1 processing

    Attributes
    ----------
        pointing_ra: float
            Pointing RA in deg
        pointing_dec: float
            Pointing DEC in deg
        target: string
            Observed source
        force_pointing: bool
            Pointing coordinates entered manualy 
            (i.e. not read from the raw FITS file)
        pointing_manual: bool
            Pointing coordinates entered manualy 
            (i.e. not read from the raw FITS file)
        input_file: string
            Input file path
        output_file: string
            Output file path
        output_file_px_charges: string
            File to store pixel charges distributions
        wobble: string
        calibration_file: string
            Calibration file used
        window_file: string
            Window file used  
        n_saturated: int
            Total number of saturated events in the file
        n_pedestal: int
            Total number of pedestal events in the file
        n_survived_pedestals
            Number of pedestal events which survived cleaning
        n_triggered_tel1: int
            Total number of TEL1 triggered events in the file
        n_triggered_tel2: int
            Total number of TEL2 triggered events in the file
        swat_event_ids_used: bool
            If true, coincident events from given night have 
            exactly the same event ID
    """

    def __init__(self):
        self.pointing_ra = float
        self.pointing_dec = float
        self.force_pointing = bool
        self.pointing_manual = bool
        self.output_file = str
        self.input_file = str
        self.output_file_px_charges = str
        self.window_file = str
        self.calibration_file = str
        self.target = str
        self.wobble = str
        self.n_pedestals = 0
        self.n_saturated = 0
        self.n_pedestals_survived = 0
        self.n_triggered_tel1 = 0
        self.n_triggered_tel2 = 0
        self.frac_rised = 0
        self.survived_charge_fraction_1 = []
        self.survived_charge_fraction_2 = []
        self.swat_event_ids_used = False

    def guess_mc(self):
        if "simtel" in self.input_file:
            return True
        else:
            return False

    def fill_target_info(self):

        target, ra_fits, dec_fits, wobble_fits = get_target(self.input_file, 
                                                            force_pointing=self.force_pointing
                                                            )
        self.target = target
        self.wobble = wobble_fits

        if (ra_fits is not None) & (dec_fits is not None):
            logging.info('Pointing info from the fits file: TARGET: ' + target + ', COORDS: ' + str(ra_fits) + ' ' + str(dec_fits) + ', WOBBLE: ' + wobble_fits)
            if self.force_pointing & (self.pointing_ra is not None) & (self.pointing_dec is not None):
                logging.info('Using pointing info from manual input anyway (forced by user).')
                self.pointing_manual=True
            else:
                self.pointing_ra = ra_fits
                self.pointing_dec = dec_fits
                self.output_file = self.output_file.split("_dl1.h5")[0] + "_" + wobble_fits + "_dl1.h5"
                self.output_file_px_charges = self.output_file_px_charges.split("_pedestal_hist.h5")[0] + "_" + wobble_fits + "_pedestal_hist.h5"
                self.pointing_manual=False

        elif self.force_pointing & (self.pointing_ra is not None) & (self.pointing_dec is not None):
            # Note that if there is no TARGET field in the fits file, it most probably means 
            # that it is a file where the shifters were tuning the trigger threshold on actual nsb conditions.
            logging.warning('No coordinates in the FITS header. Using pointing info from manual input. Are you sure that this is what you want?')
            logging.info('Pointing COORDS used: ' + str(self.pointing_ra) + ' ' + str(self.pointing_dec))
            self.pointing_manual=True

        else:
            logging.warning('No coordinates provided, exiting...')
            exit()


    def count_triggered(self, event, ismc=True):

        if ismc:
            tels = event.trigger.tels_with_trigger
            for tel in tels:
                if tel == 1:
                    self.n_triggered_tel1 += 1
                elif tel == 2:
                    self.n_triggered_tel2 += 1
        else:
            tel = event.sst1m.r0.tels_with_data[0]
            if tel == 21:
                self.n_triggered_tel1 += 1
            elif tel == 22:
                self.n_triggered_tel2 += 1

    def count_pedestals(self, event):

        tel = event.sst1m.r0.tels_with_data[0]
        event_type = event.sst1m.r0.tel[tel]._camera_event_type.value

        if event_type == 8:
            self.n_pedestals += 1

            if np.isfinite(event.dl1.tel[tel].parameters.hillas.intensity) :
                self.n_pedestals_survived += 1

    def count_survived_charge(self, event, ismc=True):

        if ismc:
            tels = event.trigger.tels_with_trigger
            for tel in tels:
                cleaning_mask = event.dl1.tel[tel].image_mask
                if tel == 1:
                    self.survived_charge_fraction_1.append(sum(event.simulation.tel[tel].true_image[cleaning_mask])/sum(event.simulation.tel[tel].true_image))
                elif tel == 2:
                    self.survived_charge_fraction_2.append(sum(event.simulation.tel[tel].true_image[cleaning_mask])/sum(event.simulation.tel[tel].true_image))
                else:
                    logging.warning('Telescope %f not recognized, survived charge fraction not logged.', tel)


    def log_result_counts(self, ismc=True):

        logging.info('Total number of TEL1 triggered events in the file: %d', 
                    self.n_triggered_tel1)
        logging.info('Total number of TEL2 triggered events in the file: %d', 
                    self.n_triggered_tel2)
        
        if not ismc:
            logging.info('Total number of saturated events in the file: %d', 
                        self.n_saturated)
            logging.info('Total number of pedestal events in the file: %d', 
                        self.n_pedestals)
            if self.n_pedestals > 0:
                logging.info('Fraction of pedestal events that survived cleaning: %f', 
                            self.n_pedestals_survived/self.n_pedestals)
            else:
                logging.info('No pedestal events found!')

