#!/bin/bash

#########################################
#  Index files of gammapy of a dataset  #
#########################################
# The DL3 files of all the nights of a dataset (submit_dataset.bash), DATASET_DIRECTORY/<night>/dl3/<analysis>/,
# in the index files of DATASET_DIRECTORY/dl3/<analysis>/ (analysis: cs1, cs2, stereo), for gammapy:
#   DataStore.from_dir("<DATASET_DIRECTORY>/dl3/stereo")

#SBATCH --job-name=dataset-dl3-index
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-01:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=8GB
#SBATCH --output=/home/%u/logs/%x-%A.out
#SBATCH --error=/home/%u/logs/%x-%A.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}/scripts/daily/common.bash"
DATASET_DIRECTORY="${DATASET_DIRECTORY:?DATASET_DIRECTORY: output directory of the dataset}"

for MAPPING in ${DATA_IRF_MODES}; do
  ANALYSIS="${MAPPING%%=*}"
  PATTERN="${DATASET_DIRECTORY}/[0-9]*/dl3/${ANALYSIS}/*.dl3.fits.gz"
  if ! compgen -G "${PATTERN}" > /dev/null; then
    echo "${ANALYSIS}: no DL3 files (${PATTERN})" >&2
    continue
  fi
  mkdir -p "${DATASET_DIRECTORY}/dl3/${ANALYSIS}"
  create_hdu_indexes --input-files "${PATTERN}" --output-dir "${DATASET_DIRECTORY}/dl3/${ANALYSIS}"
  echo "${ANALYSIS}: $(compgen -G "${PATTERN}" | wc -l) DL3 files, index in ${DATASET_DIRECTORY}/dl3/${ANALYSIS}"
done
