#!/bin/bash

#############################
#  DL1 stereo of the data   #
#############################
# Coincident events of each DL1 file of tel 1 with the DL1 files of tel 2
# (sst1mpipe_data_dl1_dl1_stereo), N_PARALLEL files at a time.

#SBATCH --job-name=data-dl1-stereo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-08:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=/home/%u/logs/%x-%j.out
#SBATCH --error=/home/%u/logs/%x-%j.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"
mkdir -p "${DATA_DIRECTORY}/dl1/stereo"

find "${DATA_DIRECTORY}/dl1/cs1" -name "*.dl1.h5" | sort | xargs -r -P "${N_PARALLEL}" -I{} \
  sst1mpipe_data_dl1_dl1_stereo --input-file-tel1 {} --input-dir-tel2 "${DATA_DIRECTORY}/dl1/cs2" \
  --config "${DATA_CONFIG}" --output-dir "${DATA_DIRECTORY}/dl1/stereo" --pattern dl1
