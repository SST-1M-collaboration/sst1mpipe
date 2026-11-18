#!/bin/bash

################################
#  DL1 -> DL2 of simulations   #
################################
# Task 3 i + j: zenith angle i of MC_ZENITHS, j = 0: point gammas, 1: diffuse gammas, 2: protons
# (test files), with the models of ${MC_DIRECTORY}/models.

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
TASK=${SLURM_ARRAY_TASK_ID:-0}
ZENITH=${ZENITHS[TASK / 3]}
NAME=${NAMES[TASK % 3]}
mkdir -p "${MC_DIRECTORY}/dl2/ze${ZENITH}"
mapfile -t MODELS < <(models_args "${MC_DIRECTORY}/models")
sst1mpipe-apply-models -i "${MC_DIRECTORY}/merged/ze${ZENITH}/${NAME}.dl1.h5" "${MODELS[@]}" \
  -o "${MC_DIRECTORY}/dl2/ze${ZENITH}/${NAME}.dl2.h5" --overwrite --no-progress
