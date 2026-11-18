#!/bin/bash

################################
#  DL1 -> DL2 of simulations   #
################################
# Task 0: point gammas, 1: diffuse gammas, 2: protons (test files), with the models of
# ${MC_DIRECTORY}/models (closest node).

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
mkdir -p "${MC_DIRECTORY}/dl2"

NAMES=(gamma_point gamma_diffuse_test proton_test)
NAME=${NAMES[${SLURM_ARRAY_TASK_ID:-0}]}
sst1mpipe-apply-models -i "${MC_DIRECTORY}/merged/${NAME}.dl1.h5" -m "${MC_DIRECTORY}/models" \
  -o "${MC_DIRECTORY}/dl2/${NAME}.dl2.h5" --overwrite --no-progress
