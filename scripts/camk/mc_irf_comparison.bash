#!/bin/bash

##############################################
#  Comparison of the IRFs of the four modes   #
##############################################
# For each zenith angle of MC_ZENITHS (task of the job array), the figures of sst1mpipe-plot-irfs with
# one series per mode (tel1, tel2, stereo), in irf/<MODELS_NAME>/<IRF_CUTS>/comparison/ze<zenith>/, with the
# performance of the SST-1M Crab paper as dashed lines of the color of the mode, where the paper has it:
# effective area at 20 deg (tel1, tel2, stereo), sensitivity, energy bias and resolution, angular resolution
# at 30 and 50 deg (tel1, stereo).

#SBATCH --job-name=mc-irf-comparison
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-01:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=/home/%u/logs/%x-%A_%a.out
#SBATCH --error=/home/%u/logs/%x-%A_%a.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"

read -ra ZENITHS <<< "${MC_ZENITHS}"
ZENITH=${ZENITHS[${SLURM_ARRAY_TASK_ID:-0}]}
IRFS="${MC_DIRECTORY}/irf/${MODELS_NAME}/${IRF_CUTS}"
DL2="${MC_DIRECTORY}/dl2/${MODELS_NAME}/ze${ZENITH}"
case "${ZENITH}" in
  20) PAPER_MODES=(tel1 tel2 stereo) ;;
  30|50) PAPER_MODES=(tel1 stereo) ;;
  *) PAPER_MODES=() ;;
esac
REFERENCES=()
for MODE in "${PAPER_MODES[@]}"; do
  REFERENCES+=(--reference "paper ${MODE}=ref:crab_paper_2025/${MODE}/ze${ZENITH}")
done
sst1mpipe-plot-irfs "${REFERENCES[@]}" \
  --irf "tel1=${IRFS}/tel1/ze${ZENITH}" --irf "tel2=${IRFS}/tel2/ze${ZENITH}" \
  --irf "stereo=${IRFS}/stereo/ze${ZENITH}" \
  --dl2 "tel1=${DL2}/tel1" --dl2 "tel2=${DL2}/tel2" --dl2 "stereo=${DL2}/stereo" \
  -o "${IRFS}/comparison/ze${ZENITH}" --title "zenith ${ZENITH} deg, ${IRF_CUTS} cuts"
