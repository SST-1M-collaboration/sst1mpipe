import glob
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
from ctapipe.core import Provenance
from ctapipe.instrument import SubarrayDescription
from ctapipe.io import DataWriter, HDF5EventSource
from ctapipe.reco import ShowerProcessor

from sst1mpipe.io import (
    get_dl1_info,
    load_dl1_pedestals,
    load_r0_pedestals,
    load_dl1_sst1m,
    load_more_dl1_tables_mono,
    write_dl1_pedestals,
    write_r0_pedestals,
)
from sst1mpipe.resources import SUBARRAY_FILE
from sst1mpipe.utils import (
    event_hillas_add_units,
    get_stereo_method,
    merge_subarrays,
    read_subarray,
    get_tel_id,
    get_telescopes,
    get_trigger_time_ns,
    hillas_prefix,
)


def find_coincidence_offset(
        tel1_file=None, tel2_files=None, outdir=None,
        config=None, save_figures=False):
    """
    Finds time offset between both telescope timestamps
    which results in the most coincident events (for
    given time window). This has to be used if the data does
    not contain White Rabbit timestamps.

    Parameters
    ----------
    tel1_file: str
        tel1 file path
    tel2_files: list of strings
    outdir: str
        Output directory path
    config: dict
    save_figures: bool

    Returns
    -------
    float

    """

    dl1_data_t1 = load_dl1_sst1m(
        tel1_file,
        tel=get_telescopes(tel1_file)[0],
        config=config,
        table='pandas',
        check_finite=True,
        quality_cuts=True
        )
    dl1_data_t2 = load_more_dl1_tables_mono(
        tel2_files,
        config=config,
        check_finite=True,
        quality_cuts=True,
        time_min=min(dl1_data_t1['local_time']),
        time_max=max(dl1_data_t1['local_time'])
        )

    timestamp_tel1 = dl1_data_t1['local_time'].to_numpy()
    timestamp_tel2 = dl1_data_t2['local_time'].to_numpy()

    # Preselection, removing clearly not coincident events to speedup the code
    mask = (timestamp_tel2 > min(timestamp_tel1) - 1.) & (timestamp_tel2 < max(timestamp_tel1) + 1.)
    timestamp_tel2_masked = timestamp_tel2[mask]
    logging.info(f'Only {len(timestamp_tel2_masked)} TEL2 events out of {len(timestamp_tel2)} kept for the offset search.')

    n_coincidences = []

    time_offsets = np.arange(
        start=config["stereo"]["SlidingWindow"]["offset_search"]["time_start"],
        stop=config["stereo"]["SlidingWindow"]["offset_search"]["time_stop"],
        step=config["stereo"]["SlidingWindow"]["offset_search"]["time_step"],
    )

    window_half_width = config["stereo"]["SlidingWindow"]["window_half_width"]

    for offset in time_offsets:

        lolim = timestamp_tel1 + offset - window_half_width
        uplim = timestamp_tel1 + offset + window_half_width

        cond_lolim = timestamp_tel2_masked >= lolim[:, np.newaxis]
        cond_uplim = timestamp_tel2_masked <= uplim[:, np.newaxis]

        mask = np.logical_and(cond_lolim, cond_uplim)

        n_coincidence = np.count_nonzero(mask)

        logging.info(f"time offset: {offset} --> {n_coincidence} events")

        n_coincidences.append([offset, n_coincidence])

    n_coincidences = np.array(n_coincidences)

    optimal_offset = n_coincidences[np.argmax(n_coincidences[:, 1]), 0]
    logging.info(f"Optimal offset found: {optimal_offset} seconds. Giving {max(n_coincidences[:, 1])} coincident events.")

    if save_figures:
        plt.figure()
        plt.plot(n_coincidences[:, 0], n_coincidences[:, 1])
        plt.grid()
        plt.xlabel('delta T [s]')
        plt.ylabel('N coincidences')
        plt.savefig(os.path.join(outdir, tel1_file.split('/')[-1].rstrip(".h5") + "_coincidences.png"))

    return optimal_offset


def make_dl1_stereo(
        dl1_file_tel1=None,
        input_dir_tel2=None,
        dl1_data_tel2=None,
        file_pattern='',
        output_path='./test.h5',
        config=None,
        time_offset=0):
    """
    Finds coincidet events in all passed tel2 DL1 files for
    given tel1 DL1 file. Resulting DL1 file contains coincident
    events only in two tabs (
    /dl1/event/telescope/parameters/{tel_021/tel_022}).
    Coincident events have the same event_id. There is also
    geometrical stereo reconstruction applied in this step
    providing extra columns to be used as features for
    stereo reconstruction (tel_impact_distance, h_max).
    This does not have to be run on MC, where the
    coincident events are found on sim_telarray level.

    Parameters
    ----------
    dl1_file_tel1: str
        tel1 file path
    input_dir_tel2: str
    dl1_data_tel2: pandas.DataFrame
    file_pattern: str
    output_path: str
    config: dict
    time_offset: float

    Returns
    -------

    """

    logging.info('Making stereo DL1 file...')
    telescope_list = []
    telescope_list.append(get_tel_id(get_telescopes(dl1_file_tel1)[0]))
    telescope_list.append(dl1_data_tel2['tel_id'][0])

    dl1_info = get_dl1_info(dl1_file_tel1)

    good_events = 0
    multiple_coincidences = 0
    mono_events = 0

    source = HDF5EventSource(dl1_file_tel1)
    current_t2_file = ""
    t2_event_id = 0

    stereo_method = get_stereo_method(config)

    if stereo_method == 'SWATEventIDs' and ~dl1_info['swat_event_ids_used']:
        logging.warning('SWATEventIDs stereo method requested, but there are no SWAT event IDs in the DL1 file. WhiteRabbitClosest method used instead.')
        stereo_method = 'WhiteRabbitClosest'

    shower_processor  = ShowerProcessor(subarray=source.subarray, config=config)

    if stereo_method == "SlidingWindow":
        window = config["stereo"]["SlidingWindow"]["window_half_width"]
    elif stereo_method == "WhiteRabbitClosest":
        window = config["stereo"]["WhiteRabbitClosest"]["max_time_diff_ns"]
    elif stereo_method == 'SWATEventIDs':
        window = 0

    source_t2 = None
    evt_t2    = None

    if stereo_method == "WhiteRabbitClosest":
        dl1_data_tel1 = load_dl1_sst1m(dl1_file_tel1, tel='tel_021', table='pandas')
        t_t1_all = get_trigger_time_ns(dl1_data_tel1)
        t_t2 = get_trigger_time_ns(dl1_data_tel2)

    # the DataWriter writes the provenance of the current activity in the output file
    Provenance().start_activity("sst1mpipe_dl1_stereo")
    with DataWriter(source,
                    output_path=output_path,
                    overwrite        = True,
                    write_dl2         = True,
                    write_dl1_parameters = True,
                    write_dl1_images         = True,
                    ) as writer:


        for _, evt in enumerate(source):

            logging.info('TEL1 event: ' + str(evt.index.event_id))

            if stereo_method == "SlidingWindow":
                t_t1 = evt.trigger.time.unix
                mask_coincidence = (dl1_data_tel2['local_time'] > t_t1 + time_offset - window) & (dl1_data_tel2['local_time'] < t_t1 + time_offset + window)
                tel2_event = dl1_data_tel2[mask_coincidence]

            elif stereo_method == "SWATEventIDs":
                ev_id1 = evt.index.event_id
                mask_coincidence = dl1_data_tel2['event_id'] == ev_id1
                tel2_event = dl1_data_tel2[mask_coincidence]

            elif stereo_method == "WhiteRabbitClosest":

                # ns trigger times of the tel1 events (trigger_time, read from the DL1 trigger table by load_dl1_sst1m)
                event_tel1_mask = dl1_data_tel1['event_id'] == evt.index.event_id

                t_t1 = t_t1_all[event_tel1_mask]
                diff = np.abs(t_t1 - t_t2)
                idx = np.argmin(diff)
                if diff[idx] < window:
                    mask_concidence = diff == diff[idx]
                    tel2_event = dl1_data_tel2[mask_concidence]
                else:
                    tel2_event = np.array([])

            tel_1 = evt.trigger.tels_with_trigger[0]
            tel1_idx = telescope_list.index(tel_1)
            tel_2 = telescope_list[tel1_idx-1]

            if tel2_event.shape[0]==1:

                if stereo_method == "SlidingWindow":
                    logging.info('Coincident event found! Time difference [s]: ' + str(abs(tel2_event['local_time'].iloc[0] - t_t1 - time_offset)))
                elif stereo_method == "WhiteRabbitClosest":
                    logging.info('Coincident event found! Time difference [ns]: ' + str(diff[idx]))
                elif stereo_method == "SWATEventIDs":
                    logging.info('Coincident event with the same SWAT Event ID found!')
                logging.info('Desired TEL ' + str(tel_2) + ' event_id: ' + str(tel2_event['event_id'].iloc[0]))

                # Looking only for coincident events to tel1 events which survived cleaning (all hillas params are calculated)
                # NOTE: Original 'hilas_reco.check_parameters(parameters=evt.dl1.tel[tel_1].parameters).all()' does not work for me
                tel1_dl1_valid = (evt.dl1.tel[tel_1].is_valid) & (evt.dl1.tel[tel_1].parameters.hillas.intensity > 0)
                tel2_dl1_valid = tel2_event[f'{hillas_prefix(tel2_event)}hillas_intensity'].iloc[0] > 0

                if tel1_dl1_valid & tel2_dl1_valid:

                    # Date and obs_id extracted from R0 filenames are stored in DL1 (DANGEROUS).
                    # obs_id should contain {DATE}{file_no}
                    datestr  = str(tel2_event['date'].iloc[0])
                    filestr  = str(tel2_event['obs_id'].iloc[0]).replace(datestr,'')

                    dl1_t2_filename = glob.glob(input_dir_tel2+f'/SST1M2*{datestr}*{filestr}*{file_pattern}*.h5')[0].split('/')[-1]
                    logging.info('Coincident events from TEL ' + str(tel_2) + ' was found in file: ' + dl1_t2_filename)
                    dl1_t2_file = os.path.join(input_dir_tel2, dl1_t2_filename)

                    evt.trigger.tels_with_trigger = np.append(evt.trigger.tels_with_trigger, tel_2)

                    # Handling a first event from the file
                    if dl1_t2_file != current_t2_file:
                        current_t2_file = dl1_t2_file
                        if source_t2 is not None:
                            source_t2.close()
                        source_t2 = iter(HDF5EventSource(current_t2_file))
                        evt_t2 = next(source_t2)

                    if tel2_event['event_id'].iloc[0] != t2_event_id:
                        t2_event_id = tel2_event['event_id'].iloc[0]

                        try:
                            while evt_t2.index.event_id != t2_event_id :
                                evt_t2 = next(source_t2)
                        except StopIteration:
                            pass
                    else:
                        logging.warning('This event from TEL ' + str(tel_2) + 'was already used! Skipping')
                        multiple_coincidences += 1
                        continue

                    evt.dl2.tel[tel_2] = evt_t2.dl2.tel[tel_2]
                    evt.dl1.tel[tel_2] = evt_t2.dl1.tel[tel_2]

                    evt.trigger.tel[tel_2].time = evt_t2.trigger.tel[tel_2].time
                    evt.monitoring.tel[tel_2].pointing = evt_t2.monitoring.tel[tel_2].pointing

                    # add units to hillas parameters (important for stereo reconstruction with shower_processor)
                    evt = event_hillas_add_units(evt)

                    shower_processor(evt)

                    writer(evt)
                    good_events += 1

                else:
                    logging.info('Coincident event in one of the two telescopes did not survived cleaning.')
            elif tel2_event.shape[0] == 0:
                logging.info('No matching event found in TEL ' + str(tel_2) + ' data. Window: ' + str(window))
                mono_events += 1
            else:
                logging.warning('More potentialy coincident events found in TEL ' + str(tel_2) + ' data. Window: ' + str(window))
                multiple_coincidences += 1

    logging.info(f"Good coincident events found: {good_events}")
    logging.info(f"Multiple coincidences: {multiple_coincidences}")
    logging.info(f"TEL {tel_1} mono events: {mono_events}")
    source.close()

    # the subarray of the tel1 DL1 file has tel1 only: the stereo file has both telescopes, at the
    # positions stored in their DL1 files (that of sst1mpipe for a telescope without file)
    tel2_files = sorted(glob.glob(os.path.join(input_dir_tel2, f'SST1M2*{file_pattern}*.h5')))
    subarrays = [read_subarray(path) for path in [dl1_file_tel1, *tel2_files[:1]]]
    subarray = merge_subarrays([*subarrays, SubarrayDescription.from_hdf(SUBARRAY_FILE)])
    subarray.select_subarray(telescope_list).to_hdf(output_path, overwrite=True)

    # Propagate tel1 pedestals
    # NOTE: We propagate only tel1 pedestals
    # Both could be stored, but I am afraid that /pedestal/tel_XX
    # table is not compatible with ctapipe readers
    logging.info('Propagating pedestals from mono tel1 DL1 to stereo DL1..')
    try:
        pedestals = load_r0_pedestals(dl1_file_tel1)
        write_r0_pedestals(output_path, pedestal_table=pedestals)
    except Exception:
        logging.warning('No pedestals found in tel1 DL1 file!')
    try:
        pedestals = load_dl1_pedestals(dl1_file_tel1)
        write_dl1_pedestals(output_path, pedestal_table=pedestals)
    except Exception:
        logging.warning('No pedestal images found in tel1 DL1 file!')
