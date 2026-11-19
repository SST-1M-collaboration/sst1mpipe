#!/bin/bash

##############################
#  R0 -> DL1 of simulations  #
##############################
# MC_FILES_PER_TASK simtel files of PARTICLE at the zenith angle ZENITH per task of the job array:
# task i, lines i * MC_FILES_PER_TASK + 1 to (i + 1) * MC_FILES_PER_TASK of
# ${MC_DIRECTORY}/lists/${PARTICLE}_ze${ZENITH}.txt (written by submit_mc.bash), or of MC_LIST (e.g.
# the files that failed).

#SBATCH --job-name=mc-r0-dl1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=2GB
#SBATCH --output=/home/%u/logs/%x-%A_%a.out
#SBATCH --error=/home/%u/logs/%x-%A_%a.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"

LIST="${MC_LIST:-${MC_DIRECTORY}/lists/${PARTICLE:?gamma_diffuse, gamma_point or proton}_ze${ZENITH:?}.txt}"
FIRST=$(( ${SLURM_ARRAY_TASK_ID:-0} * MC_FILES_PER_TASK + 1 ))
OUTPUT_DIRECTORY="${MC_DIRECTORY}/dl1/${PARTICLE}/ze${ZENITH}"
mkdir -p "${OUTPUT_DIRECTORY}"

FAILED=0
while read -r INPUT_FILE; do
  OUTPUT_FILE="${OUTPUT_DIRECTORY}/$(basename "${INPUT_FILE}" .simtel.gz).dl1.h5"
  echo "${INPUT_FILE} -> ${OUTPUT_FILE}"
  sst1mpipe-process -i "${INPUT_FILE}" -o "${OUTPUT_FILE}" --config="${MC_CONFIG}" --overwrite --no-progress \
    || { echo "Failed: ${INPUT_FILE}" >&2; FAILED=1; }
done < <(sed -n "${FIRST},$(( FIRST + MC_FILES_PER_TASK - 1 ))p" "${LIST}")
exit "${FAILED}"
