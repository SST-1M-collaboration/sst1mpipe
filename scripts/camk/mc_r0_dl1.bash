#!/bin/bash

##############################
#  R0 -> DL1 of simulations  #
##############################
# One simtel file of PARTICLE per task of the job array: line SLURM_ARRAY_TASK_ID of
# ${MC_DIRECTORY}/lists/${PARTICLE}.txt (written by submit_mc.bash).

#SBATCH --job-name=mc-r0-dl1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-08:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=2GB
#SBATCH --output=/home/%u/logs/%x-%A_%a.out
#SBATCH --error=/home/%u/logs/%x-%A_%a.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"

INPUT_FILE=$(task_line "${MC_DIRECTORY}/lists/${PARTICLE:?gamma_diffuse, gamma_point or proton}.txt")
OUTPUT_FILE="${MC_DIRECTORY}/dl1/${PARTICLE}/$(basename "${INPUT_FILE}" .simtel.gz).dl1.h5"
mkdir -p "$(dirname "${OUTPUT_FILE}")"

echo "${INPUT_FILE} -> ${OUTPUT_FILE}"
sst1mpipe-process -i "${INPUT_FILE}" -o "${OUTPUT_FILE}" --config="${MC_CONFIG}" --overwrite --no-progress
