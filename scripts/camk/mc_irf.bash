#!/bin/bash

#######################
#  IRFs (ctapipe)     #
#######################
# For each zenith angle of MC_ZENITHS and each mode, in irf/<MODELS_NAME>/<IRF_CUTS>/<mode>/ze<zenith>/:
# - tel1, tel2: one telescope alone (DL2 dl2/<MODELS_NAME>/ze<zenith>/tel1, tel2)
# - stereo: the events of the two telescopes, stereo or not (DL2 dl2/<MODELS_NAME>/ze<zenith>/stereo)
# - strict_stereo: the stereo events only (both telescopes)
# 1. ctapipe-optimize-event-selection: gammaness and direction cuts (point gammas, protons);
#    IRF_CUTS=optimized (best sensitivity) or paper (fixed: 60 % gammaness, 80 % theta)
# 2. ctapipe-compute-irf: point-like IRFs (point gammas, direction cut), and full enclosure IRFs
#    (diffuse gammas), with the benchmarks (sensitivity, resolutions)
# 3. sst1mpipe-plot-irfs: one figure per quantity in plots/
# Background: protons weighted to protons + helium (IRF_BACKGROUND, see run_irf_tool.py).
# Task 4 i + j: zenith angle i, mode j.

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
MODES=(tel1 tel2 stereo strict_stereo)
TASK=${SLURM_ARRAY_TASK_ID:-0}
ZENITH=${ZENITHS[TASK / 4]}
MODE=${MODES[TASK % 4]}

CONFIGS=()
case "${MODE}" in
  tel1|tel2) DL2="${MC_DIRECTORY}/dl2/${MODELS_NAME}/ze${ZENITH}/${MODE}" ;;
  stereo) DL2="${MC_DIRECTORY}/dl2/${MODELS_NAME}/ze${ZENITH}/stereo" ;;
  strict_stereo)
    DL2="${MC_DIRECTORY}/dl2/${MODELS_NAME}/ze${ZENITH}/stereo"
    CONFIGS+=(--config "${CONFIG_DIRECTORY}/sst1mpipe_irf_strict_stereo.yaml") ;;
esac
case "${IRF_CUTS}" in
  optimized) CUT_CONFIGS=() ;;
  paper) CUT_CONFIGS=(--config "${CONFIG_DIRECTORY}/sst1mpipe_irf_paper_cuts.yaml") ;;
  *) echo "IRF_CUTS: optimized or paper, not ${IRF_CUTS}" >&2; exit 1 ;;
esac
IRF="${MC_DIRECTORY}/irf/${MODELS_NAME}/${IRF_CUTS}/${MODE}/ze${ZENITH}"
TOOL="${SST1MPIPE_ROOT}/scripts/camk/run_irf_tool.py"
export IRF_BACKGROUND
mkdir -p "${IRF}"
echo "zenith ${ZENITH} deg, ${MODE}, ${IRF_CUTS} cuts, background ${IRF_BACKGROUND}: ${IRF}"

python "${TOOL}" ctapipe-optimize-event-selection \
  --config "${OPTIMIZE_CUTS_CONFIG}" "${CONFIGS[@]}" "${CUT_CONFIGS[@]}" \
  -g "${DL2}/gamma_point.dl2.h5" -p "${DL2}/proton_test.dl2.h5" -o "${IRF}/cuts.fits" --overwrite

python "${TOOL}" ctapipe-compute-irf --config "${COMPUTE_IRF_CONFIG}" "${CONFIGS[@]}" --cuts "${IRF}/cuts.fits" \
  --gamma-file "${DL2}/gamma_point.dl2.h5" --proton-file "${DL2}/proton_test.dl2.h5" \
  --do-background --spatial-selection-applied --DefaultFoVOffsetBins.fov_offset_n_bins=1 \
  --output "${IRF}/irf_point_like.fits.gz" --benchmark-output "${IRF}/benchmarks_point_like.fits.gz" --overwrite

python "${TOOL}" ctapipe-compute-irf --config "${COMPUTE_IRF_CONFIG}" "${CONFIGS[@]}" --cuts "${IRF}/cuts.fits" \
  --gamma-file "${DL2}/gamma_diffuse_test.dl2.h5" --proton-file "${DL2}/proton_test.dl2.h5" \
  --do-background --no-spatial-selection-applied \
  --output "${IRF}/irf_full_enclosure.fits.gz" --benchmark-output "${IRF}/benchmarks_full_enclosure.fits.gz" --overwrite

sst1mpipe-plot-irfs --irf "${MODE}=${IRF}" --dl2 "${MODE}=${DL2}" -o "${IRF}/plots" \
  --title "zenith ${ZENITH} deg, ${MODE}, ${IRF_CUTS} cuts"
