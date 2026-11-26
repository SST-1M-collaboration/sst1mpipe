#!/bin/bash

###############################################################
#  Merge the training files of several zenith angles          #
###############################################################
# Job array, one task per particle (0: gamma_diffuse, 1: proton): the training files
# merged/ze<zenith>/<particle>_train.dl1.h5 of the zenith angles of MC_TRAIN_ZENITHS (merged by
# mc_merge.bash, in this submission or before) in merged/train_<MODELS_NAME>/<particle>_train.dl1.h5.
# The obs_ids (CORSIKA runs) of each zenith angle get the offset zenith * 10000000 (offset_obs_ids.py), in
# copies of the training files: the simulations of the zenith angles have the same CORSIKA runs.

#SBATCH --job-name=mc-merge-train
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=8GB
#SBATCH --output=/home/%u/logs/%x-%A_%a.out
#SBATCH --error=/home/%u/logs/%x-%A_%a.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"

PARTICLES=(gamma_diffuse proton)
PARTICLE=${PARTICLES[${SLURM_ARRAY_TASK_ID:-0}]}
read -ra TRAIN_ZENITHS <<< "${MC_TRAIN_ZENITHS}"
OUTPUT="$(train_directory)/${PARTICLE}_train.dl1.h5"
mkdir -p "$(dirname "${OUTPUT}")"
# copies of the training files with the obs_ids of each zenith angle: obs_id + zenith * 10000000 (the
# simulations of the zenith angles have the same CORSIKA runs, which ctapipe-merge cannot merge)
FILES=()
for ZENITH in "${TRAIN_ZENITHS[@]}"; do
  FILE="${MC_DIRECTORY}/merged/ze${ZENITH}/${PARTICLE}_train.dl1.h5"
  if [[ ! -f "${FILE}" ]]; then
    echo "Missing ${FILE}: merge the zenith angle ${ZENITH} (MC_ZENITHS=\"${ZENITH}\" MC_START=merge)" >&2
    exit 1
  fi
  COPY="$(dirname "${OUTPUT}")/${PARTICLE}_train_ze${ZENITH}.tmp.dl1.h5"
  cp "${FILE}" "${COPY}"
  python "${SST1MPIPE_ROOT}/scripts/camk/offset_obs_ids.py" "${COPY}" "$(( ZENITH * 10000000 ))"
  FILES+=("${COPY}")
done
echo "${OUTPUT}: ${#FILES[@]} files"
ctapipe-merge "${FILES[@]}" --output "${OUTPUT}" --overwrite
rm -f "${FILES[@]}"
