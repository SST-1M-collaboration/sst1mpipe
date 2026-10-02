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
from collections import Counter

import astropy.units as u
import numpy as np
from ctapipe.calib import CameraCalibrator
from ctapipe.image import ImageProcessor
from ctapipe.io import DataWriter, EventSource, SimTelEventSource
from ctapipe.reco import ShowerProcessor

import sst1mpipe
from sst1mpipe.calib import (
    R0R1Calibrator,
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
from sst1mpipe.utils import (
    correct_true_image,
    energy_min_cut,
    get_swaped_modules,
    get_tel_string,
    remove_bad_pixels,
    swap_modules_r0wf,
)
from sst1mpipe.utils.monitoring_pedestals import DL1PedestalMonitor, R0PedestalMonitor, load_first_pedestals


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
                    '--max-events', '-m', type=int,
                    help='Maximum number of events read from the input file. Overrides max_events of the config file (all events if none).',
                    dest='max_events',
                    default=None
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

    outdir = args.outdir
    input_file = args.input_file
    pointing_ra = args.ra
    pointing_dec = args.dec
    force_pointing = args.force_pointing
    pixel_charges = args.pixel_charges
    reclean = args.reclean
    precise_timestamps = args.precise_timestamps

    # simtel or SST-1M zfits file, from the file content (the right EventSource is chosen by ctapipe)
    ismc = SimTelEventSource.is_compatible(input_file)

    base_name = os.path.basename(input_file)
    for suffix in (".corsika.gz.simtel.gz", ".fits.fz"):
        if base_name.endswith(suffix):
            base_name = base_name[:-len(suffix)]
    output_file = os.path.join(outdir, base_name + "_dl1.h5")
    output_logfile = os.path.join(outdir, base_name + "_r1_dl1.log")
    output_file_px_charges = os.path.join(outdir, base_name + "_pedestal_hist.h5")

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
    logging.info('Input file: %s',  input_file)
    logging.info('Output file: %s', output_file)

    # processing information written in /dl1/info
    target, wobble, pointing_manual = None, None, False
    calibration_file, window_file = None, None
    swat_event_ids_used = False
    # event counts
    n_triggered = Counter()
    n_pedestals, n_pedestals_survived, n_saturated = 0, 0, 0
    frac_rised = 0
    survived_charge_fraction = {1: [], 2: []}

    config = load_config(args.config_file, ismc=ismc)

    max_events = args.max_events if args.max_events is not None else config.get("max_events")
    if max_events is not None:
        logging.info('Maximum number of events read: %d', max_events)

    source_kwargs = {}
    if (not ismc) and force_pointing and (pointing_ra is not None) and (pointing_dec is not None):
        # pointing given by the user, used instead of the TARGET field of the file
        source_kwargs = dict(pointing_ra=pointing_ra, pointing_dec=pointing_dec)
    source = EventSource(input_url=input_file, max_events=max_events, allowed_tels=config.get("allowed_tels"), **source_kwargs)
    logging.info("Event source: %s", source.__class__.__name__)

    if source.is_simulation:
        logging.info("Tel 1 Intensity correction factor: {}".format(config['NsbCalibrator']['intensity_correction']['tel_001']))
        logging.info("Tel 2 Intensity correction factor: {}".format(config['NsbCalibrator']['intensity_correction']['tel_002']))

        if config['NsbCalibrator']['mc_correction_for_PDE']:
            used_qe = get_used_qe_simtel(source)
            logging.info("QE files used in the MC production (including the default ones): {}".format(' '.join(map(str, used_qe))))
            pde_corr_factors = get_pde_correction_factors()
            logging.info("PDE correction factors found in the calibration file mc_pde_correction_factors.json: %s", pde_corr_factors)

    else:

        logging.info("Tel 1 Intensity correction factor: {}".format(config['NsbCalibrator']['intensity_correction']['tel_021']))
        logging.info("Tel 2 Intensity correction factor: {}".format(config['NsbCalibrator']['intensity_correction']['tel_022']))

        # Target and pointing read by SST1MEventSource from the TARGET field of the Events fits header
        # (or given by the user with --force-pointing)
        target, wobble, pointing_manual = source.target, source.wobble, source.pointing_manual
        logging.info('TARGET field: %s', target)
        if (target or '').lower() in ('transition', 'dark') and not force_pointing:
            logging.info('Transition to the next wobble, or dark file, not on-source pointing direction, FILE SKIPPED.')
            exit()
        if source.pointing is None:
            logging.warning('No coordinates provided, exiting...')
            exit()
        pointing_ra, pointing_dec = source.pointing.ra.deg, source.pointing.dec.deg
        if pointing_manual:
            logging.info('Pointing COORDS used (manual input): %f %f', pointing_ra, pointing_dec)
        else:
            logging.info('Pointing info from the fits file: TARGET: %s, COORDS: %f %f, WOBBLE: %s', target, pointing_ra, pointing_dec, wobble)
            output_file = output_file.split("_dl1.h5")[0] + "_" + wobble + "_dl1.h5"
            output_file_px_charges = output_file_px_charges.split("_pedestal_hist.h5")[0] + "_" + wobble + "_pedestal_hist.h5"

        ## init the sliding windows of pedestal events and load the first pedestal events
        # r0: statistics of the ADC samples in event.mon.tel[tel].r0 (voltage drop, dead pixels)
        # dl1: statistics of the calibrated images in event.mon.tel[tel].pedestal (image cleaning)
        r0_pedestal_monitor = R0PedestalMonitor(subarray=source.subarray, config=config)
        dl1_pedestal_monitor = DL1PedestalMonitor(subarray=source.subarray, config=config)
        pedestals_in_file = load_first_pedestals(r0_pedestal_monitor, dl1_pedestal_monitor, input_file, config)

        swat_event_ids_used = source.swat_event_ids_available
        if source.swat_event_ids_available:
            logging.info('Using arrayEvtNum as event_id: input file contains SWAT array event IDs')
        else:
            logging.info('Using eventNumber as event_id: input file does not contain SWAT array event IDs')


    if reclean:
        output_file = os.path.join(outdir, output_file.split('/')[-1].rstrip(".h5") + "_recleaned.h5")
        input_file_px_charges = output_file_px_charges
        output_file_px_charges = os.path.join(outdir, output_file_px_charges.split('/')[-1].rstrip(".h5") + "_recleaned.h5")

    if source.is_simulation:
        pedestals_in_file = False

    r1_dl1_calibrator = CameraCalibrator(subarray=source.subarray, config=config)
    image_processor   = ImageProcessor(subarray=source.subarray, config=config)

    cleaner = config['ImageProcessor']['image_cleaner_type']
    # NSBImageCleaner raises the picture threshold of each pixel to pedestal_factor * std of the
    # pedestal images, taken from event.mon.tel[tel].pedestal.charge_std (see below)
    adaptive_cleaning = (cleaner == 'NSBImageCleaner') and not ismc
    pedestal_std_pe = None
    if adaptive_cleaning:
        for key, value in config['mean_charge_to_nsb_rate'].items():
            config['mean_charge_to_nsb_rate'][key] = sorted(value, key=lambda x: x['mean_charge_bin_low'], reverse=True)

        ped_mean_charge = np.ndarray(shape=[0,3])

    if reclean:
        dl1_charges = read_charge_images(input_file_px_charges)
        dl1_charges = dl1_charges[dl1_charges['n'] > config['analysis']['min_number_pedestals']]
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
        source, output_path=output_file,
        overwrite        = True,
        write_dl2         = True,
        write_dl1_parameters = True,
        write_dl1_images         = True,

    ) as writer:
        for i, event in enumerate(source):

            if not source.is_simulation:

                # NOTE: This needs to be changed in the future when event source hopefuly provides events with both telescope data
                if i == 0:
                    tel = event.trigger.tels_with_trigger[0]
                    calibrator_r0_r1 = R0R1Calibrator(subarray=source.subarray, config=config)
                    calibration_file = str(calibrator_r0_r1.calibration_file_path(tel))
                    window_corr_factors, window_file = get_window_corr_factors(
                        telescope=tel, config=config
                        )
                    tel_string = get_tel_string(tel, mc=False)
                    swaped_modules_list = get_swaped_modules(event)
                    if adaptive_cleaning:
                        dl1_pedestal_monitor.fill_monitoring(event, tel)
                        nsb_level = np.mean(event.mon.tel[tel].pedestal.charge_mean)
                        charge_to_nsb = config['mean_charge_to_nsb_rate'][tel_string]
                        for setting in charge_to_nsb:
                            min_charge = setting['mean_charge_bin_low']
                            nsb_rate = setting['nsb_rate']
                            if nsb_level >= min_charge:
                                break
                        logging.info('Average charge from the first batch of pedestal events is %f which corresponds to NSB level %s in %s', nsb_level, nsb_rate, tel_string)

                ### REAL START OF THE LOOP

                # Here we swap  wrongly connected modules
                #  swapped modules and corresponding dates
                # are stored in /data/inverted_module_list.json
                for mask_1, mask_2 in swaped_modules_list:
                    event = swap_modules_r0wf(event,mask_1, mask_2, tel=tel)

                r0_pedestal_monitor.fill_monitoring(event, tel)
                calibrator_r0_r1(event, tel)

                event_type = event.r0.tel[tel]._camera_event_type.value

                # NOTE: event.index, event.trigger and event.pointing are filled by SST1MEventSource



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

            if (not source.is_simulation) and (not reclean) and pedestals_in_file:
                # in the current setup this value is common for the whole file but keep it like this for the future
                ped_mean_charge = np.append(ped_mean_charge, [[event.index.obs_id, event.index.event_id, nsb_level]], axis=0)

            #set proper charge info according to time bins of pedestal events
            if reclean and (len(dl1_charges) > 0):
                selected_charge = dl1_charges[-1]
                for charge_entry in reversed(dl1_charges):
                    start_time = charge_entry[0]
                    if event.trigger.time >= start_time:
                        selected_charge = charge_entry
                        break
                _, _, _, sig_Qped, meanQ = selected_charge
                pedestal_std_pe = sig_Qped
                ped_mean_charge = np.append(ped_mean_charge, [[event.index.obs_id, event.index.event_id, meanQ]], axis=0)

            r1_dl1_calibrator(event) # r1->dl1a (images, peak times)

            if not source.is_simulation:

                # Integration correction of saturated pixels
                n_saturated += saturated_charge_correction(event)

                event = window_transmittance_correction(
                    event,
                    window_corr_factors=window_corr_factors,
                    telescope=tel,
                    swapped_modules=swaped_modules_list
                    )

            if adaptive_cleaning:
                # NSBImageCleaner reads the std of the pedestal images (in p.e.) from event.mon.tel[tel].pedestal
                if reclean:
                    event.mon.tel[tel].pedestal.charge_std = pedestal_std_pe
                elif pedestals_in_file:
                    # ALWAYS use adaptive cleaning - take data from the online pedestal events
                    dl1_pedestal_monitor.fill_monitoring(event, tel)
                else:
                    event.mon.tel[tel].pedestal.charge_std = None
                pedestal_std_pe = event.mon.tel[tel].pedestal.charge_std
                if pedestal_std_pe is not None:
                    picture_threshold = image_processor.clean.picture_threshold_pe.tel[tel]
                    pedestal_threshold = image_processor.clean.pedestal_factor.tel[tel] * pedestal_std_pe
                    frac_rised += np.mean(pedestal_threshold > picture_threshold)

            image_processor(event) # dl1a->dl1b (hillas parameters)

            ## Fill monitoring container with baseline info :
            if not source.is_simulation:
                if not bool(i % 100):
                    logging.info("N pixels interpolated (every 100th event): %d", calibrator_r0_r1.n_bad_pixels[tel])
                new_pedestal = False
                if event_type==8:
                    r0_pedestal_monitor(event, tel)
                    dl1_pedestal_monitor(event, tel)
                    new_pedestal = True

                elif not pedestals_in_file:

                    clenaning_mask = event.dl1.tel[tel].image_mask
                    # Arbitrary cut, just to prevent too big showers from being used
                    # We also take only every x-th event to gain some cputime
                    if (sum(clenaning_mask) < 20) and not bool(i % 10):
                        r0_pedestal_monitor(event, tel, cleaning_mask=clenaning_mask)
                        new_pedestal = True

                # writing pedestal info: ADC samples (r0) and calibrated images (dl1)
                if new_pedestal and (r0_pedestal_monitor.processed_events[tel] % 20 == 0):
                    writer._writer.write(
                        table_name='r0/monitoring/telescope/pedestal',
                        containers=[event.mon.tel[tel].r0],
                    )
                    if pedestals_in_file:
                        writer._writer.write(
                            table_name='dl1/monitoring/telescope/pedestal',
                            containers=[event.mon.tel[tel].pedestal],
                        )


            # Extraction of pixel charge distribution for MC-data tuning
            if pixel_charges:
                if not source.is_simulation:
                    event_type = event.r0.tel[tel]._camera_event_type.value
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
            n_triggered.update(event.trigger.tels_with_trigger)

            # Counting pedestal events in the file and skipping them for the output file
            if (not ismc) and event_type == 8:
                n_pedestals += 1
                n_pedestals_survived += np.isfinite(event.dl1.tel[tel].parameters.hillas.intensity)
                continue

            # Calculation of fraction of true charge which survived cleaning
            if ismc:
                for tel_id in event.trigger.tels_with_trigger:
                    true_image = event.simulation.tel[tel_id].true_image
                    if tel_id in survived_charge_fraction:
                        cleaning_mask = event.dl1.tel[tel_id].image_mask
                        survived_charge_fraction[tel_id].append(sum(true_image[cleaning_mask]) / sum(true_image))
                    else:
                        logging.warning('Telescope %d not recognized, survived charge fraction not logged.', tel_id)

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
                localtime = event.r0.tel[tel].local_camera_clock.astype(np.uint64)
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
    if (not source.is_simulation) and ((reclean and (len(dl1_charges) > 0)) or pedestals_in_file):
        write_extra_parameters(
                output_file,
                config=config, ismc=ismc, meanQ=ped_mean_charge,
                wr_timestamps=wr_timestamps
                )
    else:
        write_extra_parameters(
                output_file, config=config,
                ismc=ismc, wr_timestamps=wr_timestamps
                )

    if source.is_simulation:
        write_charge_fraction(
            output_file,
            survived_charge={
                "tel_001": survived_charge_fraction[1],
                "tel_002": survived_charge_fraction[2]
                }
            )

    # Write WR timestamps with high numerical precision
    # OBSOLETE - this function is extremely slow, there nos no reason why
    # not to extract WR timestamps in the main event loop, which makes the
    # it much faster
    #if not source.is_simulation and precise_timestamps:
    #    write_wr_timestamps(output_file,
    #                        event_source=SST1MEventSource([input_file],
    #                        max_events=max_events)
    #                        )

    # Write pointing information in the main DL1 table and in two monitoring tables
    # It is important, as we do not do it per event anymore (it was very slow)
    if not source.is_simulation:
        write_assumed_pointing(output_file, ra=pointing_ra, dec=pointing_dec, config=config)

    # Logging all event counts
    tel1_id, tel2_id = (1, 2) if ismc else (21, 22)
    logging.info('Total number of TEL1 triggered events in the file: %d', n_triggered[tel1_id])
    logging.info('Total number of TEL2 triggered events in the file: %d', n_triggered[tel2_id])
    if not ismc:
        logging.info('Total number of saturated events in the file: %d', n_saturated)
        logging.info('Total number of pedestal events in the file: %d', n_pedestals)
        if n_pedestals > 0:
            logging.info('Fraction of pedestal events that survived cleaning: %f', n_pedestals_survived / n_pedestals)
        else:
            logging.info('No pedestal events found!')

    if (not source.is_simulation) and ((reclean and (len(dl1_charges) > 0)) or pedestals_in_file):
        logging.info('Average (per event) fraction of pixels (N/1296) with raised picture threshold: %f', frac_rised/i)

    if source.is_simulation:
        # Cut on minimum mc_energy in the output file, which is needed if we want to safely combine MC from different productions
        # NOTE: This doesn't change the mc and histogram tab in the output files and this must be taken care of in performance
        # evaluation. We cannot recalculate N of simulated events at this point for each individual dl1 file, because it would
        # lead to an error of the order of 10%.
        energy_min_cut(output_file, config=config)

    # write all processing monitoring information
    write_dl1_info(output_file, dict(
        target=target, ra=pointing_ra, dec=pointing_dec, wobble=wobble, manual_coords=pointing_manual,
        calib_file=calibration_file, window_file=window_file,
        n_saturated=n_saturated, n_pedestal=n_pedestals, n_survived_pedestals=n_pedestals_survived,
        n_triggered_tel1=n_triggered[tel1_id], n_triggered_tel2=n_triggered[tel2_id],
        swat_event_ids_used=swat_event_ids_used,
    ))

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
                    write_charge_images(ped_q_map, output_file=output_file_px_charges)
            else:
                logging.warning('There are no pedestal events in the file to calculate pixel charges distributions.')
        else:
            if ((N_events_tel1 > 0) and (N_events_tel2 > 0)):
                data = np.column_stack((np.array(final_histogram_tel1), np.array(final_histogram_tel2)))
                names = ['pixel_charge_tel1', 'pixel_charge_tel2']

        if (N_events > 0) or ((N_events_tel1 > 0) and (N_events_tel2 > 0)):
            write_pixel_charges_table(data, bin_edges, names=names, output_file=output_file_px_charges)
        else:
            logging.warning('There are no pedestal events in the file to fill the pixels charge histogram.')

if __name__ == '__main__':
    main()
