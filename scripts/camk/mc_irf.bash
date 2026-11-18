#!/bin/bash

#######################
#  IRFs (ctapipe)     #
#######################
# 1. ctapipe-optimize-event-selection: gammaness and direction cuts (point gammas, protons)
# 2. ctapipe-compute-irf: point-like IRFs (point gammas, direction cut), and full enclosure IRFs
#    (diffuse gammas), with the benchmarks (sensitivity, resolutions)
# One zenith angle of MC_ZENITHS per task of the job array, in irf/ze<zenith>/.

#SBATCH --job-name=mc-irf
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
IRF="${MC_DIRECTORY}/irf/ze${ZENITH}"
mkdir -p "${IRF}"

ctapipe-optimize-event-selection --config "${OPTIMIZE_CUTS_CONFIG}" \
  -g "${DL2}/gamma_point.dl2.h5" -p "${DL2}/proton_test.dl2.h5" -o "${IRF}/cuts.fits" --overwrite

ctapipe-compute-irf --config "${COMPUTE_IRF_CONFIG}" --cuts "${IRF}/cuts.fits" \
  --gamma-file "${DL2}/gamma_point.dl2.h5" --proton-file "${DL2}/proton_test.dl2.h5" \
  --do-background --spatial-selection-applied --DefaultFoVOffsetBins.fov_offset_n_bins=1 \
  --output "${IRF}/irf_point_like.fits.gz" --benchmark-output "${IRF}/benchmarks_point_like.fits.gz" --overwrite

ctapipe-compute-irf --config "${COMPUTE_IRF_CONFIG}" --cuts "${IRF}/cuts.fits" \
  --gamma-file "${DL2}/gamma_diffuse_test.dl2.h5" --proton-file "${DL2}/proton_test.dl2.h5" \
  --do-background --no-spatial-selection-applied \
  --output "${IRF}/irf_full_enclosure.fits.gz" --benchmark-output "${IRF}/benchmarks_full_enclosure.fits.gz" --overwrite
