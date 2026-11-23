#!/bin/bash

#####################################
#  Training of the models (ctapipe) #
#####################################
# Task 0: energy regressor, 1: particle classifier, 2: disp reconstructor, in models/<MODELS_NAME>/,
# with their plots (<model>_plots): the models of the NSB level, trained on the zenith angles
# MC_TRAIN_ZENITHS (training files: train_directory)
# (feature cos_zenith).

#SBATCH --job-name=mc-train
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-12:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=/home/%u/logs/%x-%A_%a.out
#SBATCH --error=/home/%u/logs/%x-%A_%a.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"

MERGED="$(train_directory)"
MODELS="${MC_DIRECTORY}/models/${MODELS_NAME}"
mkdir -p "${MODELS}"

MODEL_NAMES=(energy_regressor particle_classifier disp_reconstructor)
MODEL=${MODEL_NAMES[${SLURM_ARRAY_TASK_ID:-0}]}
case "${MODEL}" in
  particle_classifier)
    INPUTS=(--signal "${MERGED}/gamma_diffuse_train.dl1.h5" --background "${MERGED}/proton_train.dl1.h5") ;;
  *)
    INPUTS=(-i "${MERGED}/gamma_diffuse_train.dl1.h5") ;;
esac

# shellcheck disable=SC2086  # TRAIN_OPTIONS: several options
"sst1mpipe-train-${MODEL//_/-}" "${INPUTS[@]}" -o "${MODELS}/${MODEL}.pkl" \
  --config "${CONFIG_DIRECTORY}/sst1mpipe_train_${MODEL}.yaml" --n-jobs "${SLURM_CPUS_PER_TASK:-1}" \
  --overwrite ${TRAIN_OPTIONS}
