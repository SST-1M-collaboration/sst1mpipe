#!/bin/bash

#############################
#  DL1 stereo of the data   #
#############################
# 1. ctapipe-merge: the DL1 files of the night of each telescope, dl1/merged/SST1M<tel>_<DATE>.dl1.h5
# 2. sst1mpipe-stereo: the coincident events of the two telescopes, dl1/stereo/SST1M_<DATE>.stereo.dl1.h5

#SBATCH --job-name=data-dl1-stereo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-08:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=16GB
#SBATCH --output=/home/%u/logs/%x-%j.out
#SBATCH --error=/home/%u/logs/%x-%j.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"
MERGED="${DATA_DIRECTORY}/dl1/merged"
mkdir -p "${MERGED}" "${DATA_DIRECTORY}/dl1/stereo"

for TEL in 1 2; do
  ctapipe-merge --input-dir "${DATA_DIRECTORY}/dl1/cs${TEL}" --pattern "*.dl1.h5" \
    --output "${MERGED}/SST1M${TEL}_${DATE}.dl1.h5" --overwrite --MergeTool.skip_broken_files=True
done
sst1mpipe-stereo --input "${MERGED}/SST1M1_${DATE}.dl1.h5" --input "${MERGED}/SST1M2_${DATE}.dl1.h5" \
  --output "${DATA_DIRECTORY}/dl1/stereo/SST1M_${DATE}.stereo.dl1.h5" --config "${DATA_CONFIG}" --overwrite
