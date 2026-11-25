#!/bin/bash

#####################################################
#  Comparison of the stereo combiners (camk)        #
#####################################################
# For each zenith angle of MC_ZENITHS (task of the job array), in irf/<MODELS_NAME>/<IRF_CUTS>/combiners/ze<zenith>/:
# - the figures of sst1mpipe-plot-irfs with one series per version of STEREO_VERSIONS and the SST-1M
#   Crab paper (stereo, dashed) where it has the zenith angle (submit_stereo_combiners.bash);
# - stereo_direction_<zenith>.txt/png (stereo_direction_diagnostic.py, DL2 of stereo_mars or the first
#   version): theta68 without cuts of the events with 1 or 2 telescopes, with the signs of the paper,
#   the true signs and the HillasReconstructor.

#SBATCH --job-name=mc-stereo-combiners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-02:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=16GB
#SBATCH --output=/home/%u/logs/%x-%A_%a.out
#SBATCH --error=/home/%u/logs/%x-%A_%a.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"

read -ra ZENITHS <<< "${MC_ZENITHS}"
ZENITH=${ZENITHS[${SLURM_ARRAY_TASK_ID:-0}]}
read -ra VERSIONS <<< "${STEREO_VERSIONS:-stereo_mars stereo_mars_2tel stereo_paper_1tel stereo_paper}"
IRFS="${MC_DIRECTORY}/irf/${MODELS_NAME}/${IRF_CUTS}"
DL2="${MC_DIRECTORY}/dl2/${MODELS_NAME}/ze${ZENITH}"
OUTPUT="${IRFS}/combiners/ze${ZENITH}"
mkdir -p "${OUTPUT}"

declare -A LABELS=([stereo_mars]="axes crossing, ≥1 tel" [stereo_mars_2tel]="axes crossing, 2 tel"
                   [stereo_paper_1tel]="paper signs, ≥1 tel" [stereo_paper]="paper signs, 2 tel")
ARGS=()
for VERSION in "${VERSIONS[@]}"; do
  if [[ -f "${IRFS}/${VERSION}/ze${ZENITH}/benchmarks_point_like.fits.gz" ]]; then
    ARGS+=(--irf "${LABELS[${VERSION}]:-${VERSION}}=${IRFS}/${VERSION}/ze${ZENITH}")
  else
    echo "warning: no IRFs of ${VERSION} at ${ZENITH} deg (${IRFS}/${VERSION}/ze${ZENITH})" >&2
  fi
done
case "${ZENITH}" in
  20|30|50) ARGS+=(--reference "paper ${ZENITH}°=ref:crab_paper_2025/stereo/ze${ZENITH}") ;;
esac
(( ${#ARGS[@]} > 0 )) && sst1mpipe-plot-irfs "${ARGS[@]}" -o "${OUTPUT}" --title "stereo, zenith ${ZENITH} deg, ${IRF_CUTS} cuts"

DIAGNOSTIC_VERSION=stereo_mars
[[ -f "${DL2}/${DIAGNOSTIC_VERSION}/gamma_point.dl2.h5" ]] || DIAGNOSTIC_VERSION=${VERSIONS[0]}
python "${SST1MPIPE_ROOT}/scripts/camk/stereo_direction_diagnostic.py" \
  "${DL2}/${DIAGNOSTIC_VERSION}/gamma_point.dl2.h5" "${OUTPUT}" "ze${ZENITH}"
