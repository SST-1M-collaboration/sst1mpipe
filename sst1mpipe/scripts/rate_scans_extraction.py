#!/usr/bin/env python

"""
A script to read slow control data and extract rate scans.
Based on Matthieu's code: plot_BiasCurve_fits.py
- Input is a path to a slow control directory
- Outputs are bias curves stored as txt files and png plots

Usage:

$> python rate_scans_extraction.py
--input-path /SlowControlData/
--output-dir ./
--show-figure
--only-good

"""

import argparse
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time

from sst1mpipe.io import check_outdir, load_drive_data, load_slow_data_bias_curve


def parse_args():
    parser = argparse.ArgumentParser(description="Extraction of bias curves from slow control data")

    parser.add_argument(
                        '--input-path', '-f', type=str,
                        dest='input_path',
                        help='Path to the slow control directory',
                        required=True
                        )

    parser.add_argument(
                        '--output-dir', '-o', type=str,
                        dest='outdir',
                        help='Path to store the output bias curves as txt and png.',
                        default=None
                        )
    parser.add_argument(
                    '--show-figure',
                    action='store_true',
                    help='Display figure with all rate scans in input directory.',
                    dest='show_figure'
                    )
    parser.add_argument(
                    '--only-good',
                    action='store_true',
                    help='Show and store only rate scans with associated drive files.',
                    dest='only_good'
                    )
    args = parser.parse_args()
    return args

def main():

    args = parse_args()

    aux_basepath=args.input_path
    save_path = args.outdir
    show_figure = args.show_figure
    only_good = args.only_good

    if save_path is not None:
        check_outdir(save_path)

    files = [f for f in os.listdir(aux_basepath) if os.path.isfile(os.path.join(aux_basepath, f))]

    slow_control_files = []
    drive_files = []
    for f in files:
        if 'DigicamSlowControl' in f:
            slow_control_files.append(f)

        if 'DriveSystem' in f:
            drive_files.append(f)

    if len(slow_control_files) > 0:

        font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 20}
        mpl.rc('font', **font)
        mpl.rc('legend', fontsize=15)
        plt.figure(figsize=(10, 9))

        for slow in slow_control_files:

            telescope_no=slow.split('_')[0][-1]
            file_no = telescope_no+slow.split('DigicamSlowControl'+telescope_no)[1]
            drive_file = [i for i in drive_files if file_no in i]
            slow_data = load_slow_data_bias_curve(aux_basepath+slow)

            if len(drive_file) == 1:
                drive_data = load_drive_data(aux_basepath+drive_file[0])
            else:
                print('Drive File not found!')
                drive_data = None

            save_data(slow_data, drive_data=drive_data, save_path=save_path, only_good=only_good)

        plt.yscale('log')
        plt.xlabel('Threshold [ADC]')
        plt.ylabel('Rate [Hz]')
        plt.legend()
        plt.grid()
        if save_path is not None:
            time = Time(np.array(slow_data['timestamp'][0]) / 1e3, format='unix', scale='utc')
            if save_path == ".":
                plt.savefig('bias_curve_'+str(time.to_value('iso', subfmt='date'))+'.png')
            else:
                plt.savefig(save_path+'/bias_curve_'+str(time.to_value('iso', subfmt='date'))+'.png')
        if show_figure:
            plt.show()
        plt.close()

    else:
        print('No slow control file found!')


def save_data(slow_data, drive_data=None, save_path=None, only_good=False):

    print('UTC time             Azimuth [deg]   Zenith [deg]')

    for i, val in enumerate(slow_data['biasCurvePatch7Threshold']):

        time = Time(np.array(slow_data['timestamp'][i]) / 1e3, format='unix', scale='utc')

        if drive_data is not None:
            # Now extracting average pointing position for each bias curve
            mask1 = (np.array(drive_data['TIMESTAMP']) > np.array(slow_data['timestamp'][i])-5000) * \
                    (np.array(drive_data['TIMESTAMP']) < np.array(slow_data['timestamp'][i]))

            if sum(mask1) > 0:
                avg_az = np.average(np.array(drive_data['current_position_az'])[mask1])
                avg_ze = 90-np.average(np.array(drive_data['current_position_el'])[mask1])
            else:
                avg_az = None
                avg_ze = None
        else:
            avg_az = None
            avg_ze = None
        print(time.isot, avg_az, avg_ze)

        if (avg_az is None) and only_good:
            continue

        #Now dumping to text file for comparison
        if save_path is not None:
            if save_path == ".":
                save_path = ""
            if avg_az is not None:
                outfile = save_path + 'bias_curve_' + str(slow_data['timestamp'][i] / 1e3) + '_' + str(int(avg_az)) + '_' + str(int(avg_ze)) + '.txt'
            else:
                outfile = save_path + '/bias_curve_' + str(slow_data['timestamp'][i] / 1e3) + '_no_drive.txt'
            f = open(outfile, 'w')
            f.write(f"#timestamp [s]\t{np.array(slow_data['timestamp'][i]) / 1e3}\n")
            f.write(f"#Average_azimuth\t{avg_az}\n")
            f.write(f"#Average_zenith\t{avg_ze}\n")
            f.write('#thr\ttrigger_rate\treadout_rate\tdropped_rate\n')
            for j, v in enumerate(val):
                f.write(
                    f"{v}\t"
                    f"{slow_data['biasCurveTriggerRate'][i][j]}\t"
                    f"{slow_data['biasCurveReadoutRate'][i][j]}\t"
                    f"{slow_data['biasCurveDroppedRate'][i][j]}\n"
                )
            f.close()


        plt.plot(val,  slow_data['biasCurveTriggerRate'][i], label=str(time.to_value('iso', subfmt='date_hm')))

if __name__ == '__main__':
    main()
