#!/usr/bin/env python

"""
A script to calibrate raw data (R0) or MC (R1) in DL1.
- Inputs are a single raw .fits.fz data file (containing single telescope data) 
or .simtel.gz output file of sim_telarray (may contain more telescopes).
- Output is hdf file with a table of DL1 parameters.

Usage:

$> python sst1mpipe_r0_dl1.py
--input-file SST1M1_20240304_0012.fits.fz
--output-dir ./
--config sst1mpipe_config.json
--pointing-ra 85.0
--pointing-dec 25.0
--force-pointing
--px-charges
--reclean

"""

import argparse
import logging
import os
import sys

import astropy.units as u
import numpy as np
from ctapipe.calib import CameraCalibrator
from ctapipe.io import DataWriter, EventSource
from ctapipe.reco import ShowerProcessor

import sst1mpipe
from sst1mpipe.calib import (
    Calibrator_R0_R1,
    correct_MC_for_PDE_drop,
    get_window_corr_factors,
    saturated_charge_correction,
    window_transmittance_correction,
)
from sst1mpipe.io import (
    check_outdir,
    get_pde_correction_factors,
    get_used_qe_simtel,
    load_config,
    read_charge_images,
    write_assumed_pointing,
    write_charge_fraction,
    write_charge_images,
    write_dl1_info,
    write_extra_parameters,
    write_pixel_charges_table,
)
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.utils import (
    add_event_id,
    add_pointing_to_events,
    add_trigger_time,
    correct_true_image,
    energy_min_cut,
    get_location,
    get_subarray,
    get_swap_flag,
    get_tel_string,
    image_cleaner_setup,
    remove_bad_pixels,
    swap_modules_59_88,
)
from sst1mpipe.utils.monitoring_pedestals import sliding_pedestals
from sst1mpipe.utils.monitoring_r0_dl1 import Monitoring_R0_DL1


def parse_args():

    parser = argparse.ArgumentParser(description="MC R1 or data R0 to DL1")

    # Required arguments
    parser.add_argument(
                    '--input-file', '-f', type=str,
                    dest='input_file',
                    help='Path to the simtelarray or data file',
                    required=True
                    )
    parser.add_argument('--config', '-c', action='store', type=str,
                    dest='config_file',
                    help='Path to a configuration file.',
                    required=True
                    )

    # Optional arguments
    parser.add_argument(
                    '--output-dir', '-o', type=str,
                    dest='outdir',
                    help='Path to store the output DL1 file',
                    default='./'
                    )

    parser.add_argument(
                    '--px-charges',
                    action='store_true',
                    help='Extract pixel charges for MC-data tuning and store their distribution in extra h5 file.',
                    dest='pixel_charges'
                    )

    parser.add_argument(
                    '--precise-timestamps',
                    action='store_true',
                    help='Store WR timestamps in the output DL1 table. Needs some extra processing time to go through the event source again.',
                    dest='precise_timestamps'
                    )

    parser.add_argument(
                    '--pointing-ra', '-r', type=float,
                    dest='ra',
                    help='Pointing RA (deg)',
                    default=None
                    )

    parser.add_argument(
                    '--pointing-dec', '-d', type=float,
                    dest='dec',
                    help='Pointing DEC (deg)',
                    default=None
                    )

    parser.add_argument(
                    '--force-pointing',
                    action='store_true',
                    help='Use pointing coordinates provided manualy by user even if there is a pointing info in the fits file.',
                    dest='force_pointing'
                    )

    parser.add_argument(
                    '--reclean',
                    action='store_true',
                    help='Perform cleaning based on pre-calculated charge distributions from pedestal events.',
                    dest='reclean'
                    )

    args = parser.parse_args()
    return args


def main():

    args = parse_args()
    processing_info = Monitoring_R0_DL1()

    outdir = args.outdir
    processing_info.input_file = args.input_file
    processing_info.pointing_ra = args.ra
    processing_info.pointing_dec = args.dec
    processing_info.force_pointing = args.force_pointing
    pixel_charges = args.pixel_charges
    reclean = args.reclean
    precise_timestamps = args.precise_timestamps

    ismc = processing_info.guess_mc()

    if ismc:
        processing_info.output_file = os.path.join(outdir,  processing_info.input_file.split('/')[-1].rstrip(".corsika.gz.simtel.gz") + "_dl1.h5")
        output_logfile = os.path.join(outdir, processing_info.input_file.split('/')[-1].rstrip(".corsika.gz.simtel.gz") + "_r1_dl1.log")
        processing_info.output_file_px_charges = os.path.join(outdir, processing_info.input_file.split('/')[-1].rstrip(".corsika.gz.simtel.gz") + "_pedestal_hist.h5")
    else:
        processing_info.output_file = os.path.join(outdir,  processing_info.input_file.split('/')[-1].rstrip(".fits.fz") + "_dl1.h5")
        output_logfile = os.path.join(outdir,  processing_info.input_file.split('/')[-1].rstrip(".fits.fz") + "_r1_dl1.log")
        processing_info.output_file_px_charges = os.path.join(outdir,  processing_info.input_file.split('/')[-1].rstrip(".fits.fz") + "_pedestal_hist.h5")

    check_outdir(outdir)

    if reclean:
        output_logfile = os.path.join(outdir, output_logfile.split('/')[-1].rstrip(".log") + "_recleaned.log")

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers= [
            logging.FileHandler(output_logfile, 'w+'),
            logging.StreamHandler(stream=sys.stdout)
            ]
    )

    logging.info('sst1mpipe version: %s', sst1mpipe.__version__)
    logging.info('Input file: %s',  processing_info.input_file)
    logging.info('Output file: %s', processing_info.output_file)

    max_events = None

    config = load_config(args.config_file, ismc=ismc)

    if ismc:
        source = EventSource(processing_info.input_file, max_events=max_events, allowed_tels=config["allowed_tels"])

        logging.info("Tel 1 Intensity correction factor: {}".format(config['NsbCalibrator']['intensity_correction']['tel_001']))
        logging.info("Tel 2 Intensity correction factor: {}".format(config['NsbCalibrator']['intensity_correction']['tel_002']))

        if config['NsbCalibrator']['mc_correction_for_PDE']:
            used_qe = get_used_qe_simtel(source)
            logging.info("QE files used in the MC production (including the default ones): {}".format(' '.join(map(str, used_qe))))
            pde_corr_factors = get_pde_correction_factors()
            logging.info("PDE correction factors found in the calibration file mc_pde_correction_factors.json: %s", pde_corr_factors)

    else:
        source = SST1MEventSource([processing_info.input_file], max_events=max_events)
        source._subarray = get_subarray()

        logging.info("Tel 1 Intensity correction factor: {}".format(config['NsbCalibrator']['intensity_correction']['tel_021']))
        logging.info("Tel 2 Intensity correction factor: {}".format(config['NsbCalibrator']['intensity_correction']['tel_022']))

        ## init pedestal_info and loading first pedestal events in pedestal_info
        pedestal_info = sliding_pedestals(input_file=processing_info.input_file, config=config)

        # Reading target name and assumed pointing ra,dec from the target field
        # of the Events fits header
        processing_info.fill_target_info()

        processing_info.swat_event_ids_used = source.swat_event_ids_available
        if source.swat_event_ids_available:
            logging.info('Using arrayEvtNum as event_id: input file contains SWAT array event IDs')
        else:
            logging.info('Using eventNumber as event_id: input file does not contain SWAT array event IDs')
    

    if reclean:
        processing_info.output_file = os.path.join(outdir, processing_info.output_file.split('/')[-1].rstrip(".h5") + "_recleaned.h5")
        input_file_px_charges = processing_info.output_file_px_charges
        processing_info.output_file_px_charges = os.path.join(outdir, processing_info.output_file_px_charges.split('/')[-1].rstrip(".h5") + "_recleaned.h5")

    r1_dl1_calibrator = CameraCalibrator(subarray=source.subarray, config=config)
    image_processor   = image_cleaner_setup(subarray=source.subarray, config=config, ismc=ismc)

    cleaner = config['ImageProcessor']['image_cleaner_type']
    if (cleaner == 'ImageCleanerSST') and not ismc:
        #to be sure about the order, sort telescope_defaults according to increasing nsb level
        for key, value in config['ImageProcessor'][cleaner]['telescope_defaults'].items():
            config['ImageProcessor'][cleaner]['telescope_defaults'][key] = sorted(value, key=lambda x: x['min_nsb_level'], reverse=True)

        for key, value in config['mean_charge_to_nsb_rate'].items():
            config['mean_charge_to_nsb_rate'][key] = sorted(value, key=lambda x: x['mean_charge_bin_low'], reverse=True)

        image_processor.clean.nsb_level = np.mean(pedestal_info.get_img_charge_mean())
        image_processor.clean.config = config['ImageProcessor'][cleaner]
        ped_mean_charge = np.ndarray(shape=[0,3])

    if reclean:
        dl1_charges = read_charge_images(input_file_px_charges)
        dl1_charges = dl1_charges[dl1_charges['n'] > config['ImageProcessor'][cleaner]['min_number_pedestals']]
        dl1_charges['mean_charge'] = np.average(dl1_charges['average_q'], axis=1)

    shower_processor  = ShowerProcessor(subarray=source.subarray, config=config)

    if pixel_charges:
        BINS = 1000
        N_events= 0
        N_events_tel1= 0
        N_events_tel2= 0
        if not source.is_simulation:
            final_histogram = np.zeros(BINS)
            #NOTE: To store images of stdevs. and average charges for pedestal events
            ped_q_map = []
            ped_q_sum = 0
            ped_q2_sum = 0
            ped_n = 0
            ped_time_start = None
            ped_time_window = config["analysis"]["ped_time_window"]*u.s
        else:
            final_histogram_tel1 = np.zeros(BINS)
            final_histogram_tel2 = np.zeros(BINS)

    if not source.is_simulation and precise_timestamps:
        full_seconds = []
        fractional_seconds = []

    with DataWriter(
        source, output_path=processing_info.output_file, 
        overwrite        = True, 
        write_showers    = True,
        write_parameters = True,
        write_images     = True,

    ) as writer:
        for i, event in enumerate(source):

            if not source.is_simulation:

                # NOTE: This needs to be changed in the future when event source hopefuly provides events with both telescope data
                if i == 0:
                    tel = event.sst1m.r0.tels_with_data[0]
                    calibrator_r0_r1 = Calibrator_R0_R1(config=config, telescope=tel)
                    processing_info.calibration_file = calibrator_r0_r1.calibration_file
                    window_corr_factors, processing_info.window_file = get_window_corr_factors(
                        telescope=tel, config=config
                        )
                    tel_string = get_tel_string(tel, mc=False)
                    location = get_location(config=config, tel=tel_string)
                    swap_modules = get_swap_flag(event)
                    if (cleaner == 'ImageCleanerSST'):
                        charge_to_nsb = config['mean_charge_to_nsb_rate'][tel_string]
                        for setting in charge_to_nsb:
                            min_charge = setting['mean_charge_bin_low']
                            nsb_rate = setting['nsb_rate']
                            if image_processor.clean.nsb_level >= min_charge:
                                break
                        logging.info('Average charge from the first batch of pedestal events is %f which corresponds to NSB level %s in %s', image_processor.clean.nsb_level, nsb_rate, tel_string)

                event.trigger.tels_with_trigger = [tel]

                event = calibrator_r0_r1.calibrate(event, pedestal_info=pedestal_info)
                # print(calibrator_r0_r1.pixels_removed) # can be monitored

                event.r1.tel[tel].selected_gain_channel = np.zeros(source.subarray.tel[tel].camera.readout.n_pixels,dtype='int8')
                event_type = event.sst1m.r0.tel[tel]._camera_event_type.value

                # Fill trigger container properly
                # SHOULD BE REMOVED as soon as event source can handle this
                event = add_trigger_time(event, telescope=tel)

                # Add assumed pointing (this should be part of Event Source in the future)
                # This stores the pointing information in the right containters. If done this way, pointing information is automaticaly propagated in
                # the output DL1 file, in /dl1/monitoring/subarray/pointing and /dl1/monitoring/telescope/pointing/TEL
                # This takes absurdly long time - revisit in the future, checking if we can make it faster
                event = add_pointing_to_events(
                                                event, 
                                                ra=processing_info.pointing_ra, 
                                                dec=processing_info.pointing_dec, 
                                                telescope=tel, 
                                                location=location
                                                )

                # Adding event_id and obs_id in event.index
                # event_id is in event.sst1m.r0.event_id, but obs_id must be made up
                # SHOULD BE REMOVED as soon as event source can handle this
                event = add_event_id(event, filename=processing_info.output_file, event_number=i)

                # Here we swap two wrongly connected modules in tel2 after 
                # camera refurbishment in 2023. Only R1 waveforms are swapped.
                # The modules were physicaly reconnected in July 2024, so no swapping
                # after this date (authomatic)
                event = swap_modules_59_88(event, tel=tel, swap_flag=swap_modules)

            # For an unknown reason, event.simulation.tel[tel].true_image is sometime None, which kills the rest of the script
            # and simulation histogram is not saved. Here we repace it with an array of zeros.
            if source.is_simulation:
                event = correct_true_image(event)
                # Now include PDE correction based on the PDE drop set in MC
                if config['NsbCalibrator']['mc_correction_for_PDE']:
                    event = correct_MC_for_PDE_drop(event, 
                        simtel_config_qe=used_qe,
                        pde_corr_factors=pde_corr_factors
                        )

            # This function flags the bad pixel according to the cfg file, and just for sure also kills the waveforms.
            # Charges in these pixels are then interpolated using method set in cfg: invalid_pixel_handler_type
            # Default is NeighborAverage, but can be turned off with 'null'
            event = remove_bad_pixels(event, config=config)
            
            if (not source.is_simulation) and (not reclean) and pedestal_info.pedestals_in_file:
                # ALWAYS use adaptive cleaning - take data from online pedestal_info
                image_processor.clean.average_charge = pedestal_info.get_img_charge_mean()
                image_processor.clean.stdev_charge = pedestal_info.get_img_charge_std()
                # in the current setup this value is common for the whole file but keep it like this for the future
                ped_mean_charge = np.append(ped_mean_charge, [[event.index.obs_id, event.index.event_id, image_processor.clean.nsb_level]], axis=0)
            
            #set proper charge info according to time bins of pedestal events
            if reclean and (len(dl1_charges) > 0):
                for [start_time, n, Qped, sig_Qped, meanQ] in reversed(dl1_charges):
                    if event.trigger.time >= start_time:
                        break
                image_processor.clean.average_charge = Qped
                image_processor.clean.stdev_charge = sig_Qped
                image_processor.clean.nsb_level = meanQ
                image_processor.clean.config = config['ImageProcessor'][cleaner]
                ped_mean_charge = np.append(ped_mean_charge, [[event.index.obs_id, event.index.event_id, meanQ]], axis=0)

            r1_dl1_calibrator(event) # r1->dl1a (images, peak times)                

            if not source.is_simulation:

                # Integration correction of saturated pixels
                event = saturated_charge_correction(event, processing_info=processing_info)

                event = window_transmittance_correction(
                    event, 
                    window_corr_factors=window_corr_factors, 
                    telescope=tel,
                    swap_flag=swap_modules
                    )

            image_processor(event) # dl1a->dl1b (hillas parameters)

            if (cleaner == 'ImageCleanerSST') and not ismc:
                processing_info.frac_rised += image_processor.clean.frac_rised

            ## Fill monitoring container with baseline info :
            if not source.is_simulation:
                if not bool(i % 100) and config["telescope_calibration"]["bad_calib_px_interpolation"]:
                    logging.info("N pixels interpolated (every 100th event): %d", calibrator_r0_r1.pixels_removed)
                if event_type==8:
                    pedestal_info.add_ped_evt(event)

                    # writing pedestal info in dl1
                    if ( (pedestal_info.processed_pedestals !=0) and \
                         (pedestal_info.processed_pedestals%20 == 0)):

                        pedestal_info.fill_mon_container(event)
                        writer._writer.write(
                            table_name='dl1/monitoring/telescope/pedestal',
                            containers=[event.mon.tel[tel].pedestal],
                        )

                elif not pedestal_info.pedestals_in_file:

                    clenaning_mask = event.dl1.tel[tel].image_mask
                    # Arbitrary cut, just to prevent too big showers from being used
                    # We also take only every x-th event to gain some cputime
                    if (sum(clenaning_mask) < 20) and not bool(i % 10):
                        pedestal_info.add_ped_evt(event, cleaning_mask=clenaning_mask, store_image=False)

                    # writing pedestal info in dl1
                    if ( (pedestal_info.processed_pedestals !=0) and \
                         (pedestal_info.processed_pedestals%20 == 0)):

                        pedestal_info.fill_mon_container(event)
                        writer._writer.write(
                            table_name='dl1/monitoring/telescope/pedestal',
                            containers=[event.mon.tel[tel].pedestal],
                        )


            # Extraction of pixel charge distribution for MC-data tuning
            if pixel_charges:
                if not source.is_simulation:
                    event_type = event.sst1m.r0.tel[tel]._camera_event_type.value
                    if ped_time_start is None:
                        ped_time_start = event.trigger.time
                    if event_type == 8:
                        image = event.dl1.tel[tel].image
                        hist, bin_edges = np.histogram(image, range=(-10, 40), bins=BINS, density=False)
                        final_histogram += hist
                        N_events += 1
                        time_diff = event.trigger.time - ped_time_start
                        if (time_diff > ped_time_window) and (ped_n > 0):
                            mean = ped_q_sum/ped_n
                            stdev = np.sqrt(ped_q2_sum/ped_n-mean*mean)
                            ped_q_map.append([ped_time_start, ped_n, mean, stdev])
                            ped_n = 0
                            ped_q_sum = 0
                            ped_q2_sum = 0
                            ped_time_start = event.trigger.time
                        ped_q_sum += image
                        ped_q2_sum += image*image
                        ped_n += 1

                else:
                    for tel in event.trigger.tels_with_trigger:
                        # We need to get rid of shower pixels
                        noise_mask = ~np.array(event.simulation.tel[tel].true_image, dtype=bool)
                        image = event.dl1.tel[tel].image[noise_mask]
                        hist, bin_edges = np.histogram(image, range=(-10, 40), bins=BINS, density=False)
                        if tel == 1:
                            final_histogram_tel1 += hist
                            N_events_tel1 += 1
                        elif tel == 2:
                            final_histogram_tel2 += hist
                            N_events_tel2 += 1

            # We would like to store in DL1 also some additional parameters needed for disp reconstruction and few more additional features
            # It cannot be done at this level, because: AttributeError: 'CameraHillasParametersContainer' object has no attribute 'disp'

            shower_processor(event) # dl1b->dl2 (reconstruction of stereo parameters, also energy/direction/classification in the future versions of ctapipe)

            # Counting all triggered events
            processing_info.count_triggered(event, ismc=source.is_simulation)

            # Counting pedestal events in the file and skipping them for the output file
            if not ismc:
                processing_info.count_pedestals(event)

            # skip rest of the script for pedestal events
            if not ismc:
                if event_type == 8:
                    continue

            # Calculation of fraction of true charge which survived cleaning
            processing_info.count_survived_charge(event, ismc=source.is_simulation)
            
            ## Correct (or not) the Voltage drop effect : Global correction on the intensity
            ## apply (or not) some absolute correction on the intensity
            
            if not source.is_simulation:
                I0 = event.dl1.tel[tel].parameters.hillas.intensity
                # VN: to be consistent during the cleaning I had to move all gain drop corrections into calibration
                I_corr = I0*config['NsbCalibrator']["intensity_correction"][tel_string]
                event.dl1.tel[tel].parameters.hillas.intensity = I_corr
            else:
                for tel in event.trigger.tels_with_trigger:
                    tel_string = get_tel_string(tel, mc=True)
                    I0 = event.dl1.tel[tel].parameters.hillas.intensity
                    I_corr = I0*config['NsbCalibrator']["intensity_correction"][tel_string]
                    event.dl1.tel[tel].parameters.hillas.intensity = I_corr

            writer(event)

            # Extracting WR timestamps with high numerical precision
            if not source.is_simulation and precise_timestamps:
                localtime = event.sst1m.r0.tel[tel].local_camera_clock.astype(np.uint64)
                S_TO_NS = np.uint64(1e9)
                full_seconds.append(localtime // S_TO_NS)
                fractional_seconds.append((localtime % S_TO_NS) / S_TO_NS)


        if max_events is None and source.is_simulation:
            writer.write_simulation_histograms(source)

    if not source.is_simulation and precise_timestamps:
        wr_timestamps = np.column_stack((full_seconds, fractional_seconds))
    else: 
        wr_timestamps=None

    # Write additional params in the DL1 file
    # - these are not defined in the ctapipe containers, but are necessary for (mono) reconstruction
    # - DISP parameters are calculated and stored
    # - some more parameters are extracted from other tables in the file and added to the parameters table for convenience
    # NOTE: unfortunately using this, units in all columns have to be dropped, because otherwise ctapipe merging tool fails.
    # I didn't find a solution, this should definitely be revisited!
    if (not source.is_simulation) and ((reclean and (len(dl1_charges) > 0)) or pedestal_info.pedestals_in_file):
        write_extra_parameters(
                processing_info.output_file, 
                config=config, ismc=ismc, meanQ=ped_mean_charge, 
                wr_timestamps=wr_timestamps
                )
    else:
        write_extra_parameters(
                processing_info.output_file, config=config, 
                ismc=ismc, wr_timestamps=wr_timestamps
                )

    if source.is_simulation:
        write_charge_fraction(
            processing_info.output_file, 
            survived_charge={
                "tel_001": processing_info.survived_charge_fraction_1, 
                "tel_002": processing_info.survived_charge_fraction_2
                }
            )

    # Write WR timestamps with high numerical precision
    # OBSOLETE - this function is extremely slow, there nos no reason why
    # not to extract WR timestamps in the main event loop, which makes the
    # it much faster
    #if not source.is_simulation and precise_timestamps:
    #    write_wr_timestamps(processing_info.output_file, 
    #                        event_source=SST1MEventSource([processing_info.input_file], 
    #                        max_events=max_events)
    #                        )

    # Write pointing information in the main DL1 table and in two monitoring tables
    # It is important, as we do not do it per event anymore (it was very slow)
    if not source.is_simulation:
        write_assumed_pointing(processing_info, config=config)

    # Logging all event counts
    processing_info.log_result_counts(ismc=source.is_simulation)

    if (not source.is_simulation) and ((reclean and (len(dl1_charges) > 0)) or pedestal_info.pedestals_in_file):
        logging.info('Average (per event) fraction of pixels (N/1296) with raised picture threshold: %f', processing_info.frac_rised/i)
        # to dump how many times particular pixels were raised
        #image_processor.clean.dump()

    if source.is_simulation:
        # Cut on minimum mc_energy in the output file, which is needed if we want to safely combine MC from different productions
        # NOTE: This doesn't change the mc and histogram tab in the output files and this must be taken care of in performance
        # evaluation. We cannot recalculate N of simulated events at this point for each individual dl1 file, because it would 
        # lead to an error of the order of 10%.
        energy_min_cut(processing_info.output_file, config=config)

    # write all processing monitoring information
    write_dl1_info(processing_info)

    # We write calibration configuration in the output file
    # NOTE: If one use the ctapipe merging tool this table is missing in the merged DL1 file!
    # TODO: Broken after implementation telescope dependent tailcuts, but not supper important
    # write_r1_dl1_cfg(output_file, config=config)

    # Save pixel charges histograms and maps in output file
    if pixel_charges:

        if not source.is_simulation:

            if ped_n > 0:
                mean = ped_q_sum/ped_n
                stdev = np.sqrt(ped_q2_sum/ped_n-mean*mean)
                ped_q_map.append([ped_time_start, ped_n, mean, stdev])

                if (N_events > 0):
                    data = np.array(final_histogram)[..., np.newaxis]
                    names = ['pixel_charge']
                    write_charge_images(ped_q_map, output_file=processing_info.output_file_px_charges)
            else:
                logging.warning('There are no pedestal events in the file to calculate pixel charges distributions.')
        else:
            if ((N_events_tel1 > 0) and (N_events_tel2 > 0)):
                data = np.column_stack((np.array(final_histogram_tel1), np.array(final_histogram_tel2)))
                names = ['pixel_charge_tel1', 'pixel_charge_tel2']
        
        if (N_events > 0) or ((N_events_tel1 > 0) and (N_events_tel2 > 0)):
            write_pixel_charges_table(data, bin_edges, names=names, output_file=processing_info.output_file_px_charges)
        else:
            logging.warning('There are no pedestal events in the file to fill the pixels charge histogram.')

if __name__ == '__main__':
    main()
