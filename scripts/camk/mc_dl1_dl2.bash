#!/bin/bash

################################
#  DL1 -> DL2 of simulations   #
################################
# Test files (point gammas, diffuse gammas, protons) of each zenith angle of MC_ZENITHS, with the
# models of ${MC_DIRECTORY}/models, in three versions: dl2/ze<zenith>/<version>/
# - stereo: both telescopes (stereo combination, and the events of one telescope)
# - tel1, tel2: one telescope alone (its mono reconstruction)
# Task 9 i + 3 j + k: zenith angle i, version j, file k.

#SBATCH --job-name=mc-dl1-dl2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-08:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=8GB
#SBATCH --output=/home/%u/logs/%x-%A_%a.out
#SBATCH --error=/home/%u/logs/%x-%A_%a.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"

read -ra ZENITHS <<< "${MC_ZENITHS}"
NAMES=(gamma_point gamma_diffuse_test proton_test)
VERSIONS=(stereo tel1 tel2)
declare -A TELESCOPES=([stereo]="[]" [tel1]="[1]" [tel2]="[2]")
TASK=${SLURM_ARRAY_TASK_ID:-0}
ZENITH=${ZENITHS[TASK / 9]}
VERSION=${VERSIONS[TASK % 9 / 3]}
NAME=${NAMES[TASK % 3]}
OUTPUT="${MC_DIRECTORY}/dl2/ze${ZENITH}/${VERSION}/${NAME}.dl2.h5"
mkdir -p "$(dirname "${OUTPUT}")"
mapfile -t MODELS < <(models_args "${MC_DIRECTORY}/models")
sst1mpipe-apply-models -i "${MC_DIRECTORY}/merged/ze${ZENITH}/${NAME}.dl1.h5" "${MODELS[@]}" -o "${OUTPUT}" \
  --ApplyModels.telescopes="${TELESCOPES[${VERSION}]}" --overwrite --no-progress
