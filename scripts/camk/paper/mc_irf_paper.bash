#!/bin/bash

################################################################
#  IRFs with the cuts of the SST-1M Crab paper (comparison)    #
################################################################
# As mc_irf.bash, for the comparison with arXiv:2506.01733: fixed cuts (60 % gamma efficiency,
# 80 % theta containment), stereo events only, protons weighted to protons + helium (DAMPE).
# One zenith angle of MC_ZENITHS per task of the job array, in irf_paper/ze<zenith>/.

#SBATCH --job-name=mc-irf-paper
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-04:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=8GB
#SBATCH --output=/home/%u/logs/%x-%A_%a.out
#SBATCH --error=/home/%u/logs/%x-%A_%a.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"

read -ra ZENITHS <<< "${MC_ZENITHS}"
ZENITH=${ZENITHS[${SLURM_ARRAY_TASK_ID:-0}]}
DL2="${MC_DIRECTORY}/dl2/ze${ZENITH}"
IRF="${MC_DIRECTORY}/irf_paper/ze${ZENITH}"
PAPER="${SST1MPIPE_ROOT}/scripts/camk/paper"
mkdir -p "${IRF}"

python "${PAPER}/run_p_he.py" ctapipe-optimize-event-selection \
  --config "${OPTIMIZE_CUTS_CONFIG}" --config "${PAPER}/irf_paper_cuts.yaml" \
  -g "${DL2}/gamma_point.dl2.h5" -p "${DL2}/proton_test.dl2.h5" -o "${IRF}/cuts.fits" --overwrite

python "${PAPER}/run_p_he.py" ctapipe-compute-irf \
  --config "${COMPUTE_IRF_CONFIG}" --config "${PAPER}/irf_paper_cuts.yaml" --cuts "${IRF}/cuts.fits" \
  --gamma-file "${DL2}/gamma_point.dl2.h5" --proton-file "${DL2}/proton_test.dl2.h5" \
  --do-background --spatial-selection-applied --DefaultFoVOffsetBins.fov_offset_n_bins=1 \
  --output "${IRF}/irf_point_like.fits.gz" --benchmark-output "${IRF}/benchmarks_point_like.fits.gz" --overwrite
