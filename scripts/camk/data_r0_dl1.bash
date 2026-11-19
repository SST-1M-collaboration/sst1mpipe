#!/bin/bash

#############################
#  R0 -> DL1 of the data    #
#############################
# One raw file of telescope TEL (1 or 2) per task of the job array: line SLURM_ARRAY_TASK_ID of
# ${DATA_DIRECTORY}/lists/cs${TEL}.txt (written by submit_data.bash). sst1mpipe-monitoring, then
# sst1mpipe-process with the monitoring file (the wobble is added to the DL1 file name).

#SBATCH --job-name=data-r0-dl1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-04:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=/home/%u/logs/%x-%A_%a.out
#SBATCH --error=/home/%u/logs/%x-%A_%a.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"

INPUT_FILE=$(task_line "${DATA_DIRECTORY}/lists/cs${TEL:?1 or 2}.txt")
RUN=$(basename "${INPUT_FILE}" .fits.fz)
MONITORING_FILE="${DATA_DIRECTORY}/monitoring/cs${TEL}/${RUN}.monitoring.h5"
OUTPUT_FILE="${DATA_DIRECTORY}/dl1/cs${TEL}/${RUN}.dl1.h5"
mkdir -p "$(dirname "${MONITORING_FILE}")" "$(dirname "${OUTPUT_FILE}")"

echo "${INPUT_FILE} -> ${OUTPUT_FILE}"
sst1mpipe-monitoring -i "${INPUT_FILE}" -o "${MONITORING_FILE}" --config "${MONITORING_CONFIG}" --overwrite
sst1mpipe-process -i "${INPUT_FILE}" -o "${OUTPUT_FILE}" --config "${DATA_CONFIG}" \
  --monitoring-source SST1MMonitoringSource --SST1MMonitoringSource.input_files "${MONITORING_FILE}" \
  --overwrite --no-progress
