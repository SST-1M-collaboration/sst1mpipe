#!/bin/bash

#######################
#  IRFs (ctapipe)     #
#######################
# For each zenith angle of MC_ZENITHS and each mode, in irf/<MODELS_NAME>/<IRF_CUTS>/<mode>/ze<zenith>/:
# - tel1, tel2: one telescope alone (DL2 dl2/<MODELS_NAME>/ze<zenith>/tel1, tel2)
# - stereo: the events with both telescopes (stereo models, DL2 dl2/<MODELS_NAME>/ze<zenith>/stereo)
# 1. sst1mpipe-optimize-event-selection: gammaness and direction cuts (point gammas, protons);
#    IRF_CUTS=optimized (best sensitivity) or paper (60 % gammaness per energy bin, one theta cut of 80 %)
# 2. sst1mpipe-compute-irf: point-like IRFs (point gammas, direction cut), and full enclosure IRFs
#    (diffuse gammas), with the benchmarks (sensitivity, resolutions)
# 3. sst1mpipe-plot-irfs: one figure per quantity in plots/
# Spectra of the SST-1M Crab paper: Crab Nebula of MAGIC, protons weighted to protons + helium
# (IRF_BACKGROUND=p_he, or protons; sst1mpipe.performance.irf_tools).
# Task 3 i + j: zenith angle i, mode j.

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
MODES=(tel1 tel2 stereo)
TASK=${SLURM_ARRAY_TASK_ID:-0}
ZENITH=${ZENITHS[TASK / 3]}
MODE=${MODES[TASK % 3]}

CONFIGS=()
DL2="${MC_DIRECTORY}/dl2/${MODELS_NAME}/ze${ZENITH}/${MODE}"
case "${IRF_CUTS}" in
  optimized) CUT_CONFIGS=() ;;
  paper) CUT_CONFIGS=(--config "${CONFIG_DIRECTORY}/sst1mpipe_irf_paper_cuts.yaml") ;;
  *) echo "IRF_CUTS: optimized or paper, not ${IRF_CUTS}" >&2; exit 1 ;;
esac
IRF="${MC_DIRECTORY}/irf/${MODELS_NAME}/${IRF_CUTS}/${MODE}/ze${ZENITH}"
mkdir -p "${IRF}"
echo "zenith ${ZENITH} deg, ${MODE}, ${IRF_CUTS} cuts, background ${IRF_BACKGROUND}: ${IRF}"

sst1mpipe-optimize-event-selection --background "${IRF_BACKGROUND}" \
  --config "${OPTIMIZE_CUTS_CONFIG}" "${CONFIGS[@]}" "${CUT_CONFIGS[@]}" \
  -g "${DL2}/gamma_point.dl2.h5" -p "${DL2}/proton_test.dl2.h5" -o "${IRF}/cuts.fits" --overwrite

sst1mpipe-compute-irf --background "${IRF_BACKGROUND}" --config "${COMPUTE_IRF_CONFIG}" "${CONFIGS[@]}" --cuts "${IRF}/cuts.fits" \
  --gamma-file "${DL2}/gamma_point.dl2.h5" --proton-file "${DL2}/proton_test.dl2.h5" \
  --do-background --spatial-selection-applied --DefaultFoVOffsetBins.fov_offset_n_bins=1 \
  --output "${IRF}/irf_point_like.fits.gz" --benchmark-output "${IRF}/benchmarks_point_like.fits.gz" --overwrite

sst1mpipe-compute-irf --background "${IRF_BACKGROUND}" --config "${COMPUTE_IRF_CONFIG}" "${CONFIGS[@]}" --cuts "${IRF}/cuts.fits" \
  --gamma-file "${DL2}/gamma_diffuse_test.dl2.h5" --proton-file "${DL2}/proton_test.dl2.h5" \
  --do-background --no-spatial-selection-applied \
  --output "${IRF}/irf_full_enclosure.fits.gz" --benchmark-output "${IRF}/benchmarks_full_enclosure.fits.gz" --overwrite

sst1mpipe-plot-irfs --irf "${MODE}=${IRF}" --dl2 "${MODE}=${DL2}" -o "${IRF}/plots" \
  --title "zenith ${ZENITH} deg, ${MODE}, ${IRF_CUTS} cuts"
