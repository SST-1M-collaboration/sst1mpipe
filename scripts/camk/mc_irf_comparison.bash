#!/bin/bash

##############################################
#  Comparison of the IRFs of the four modes   #
##############################################
# For each zenith angle of MC_ZENITHS (task of the job array), the figures of sst1mpipe-plot-irfs with
# one series per mode (tel1, tel2, stereo, strict_stereo), in irf/<IRF_CUTS>/comparison/ze<zenith>/.

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
IRFS="${MC_DIRECTORY}/irf/${IRF_CUTS}"
DL2="${MC_DIRECTORY}/dl2/ze${ZENITH}"
sst1mpipe-plot-irfs \
  --irf "tel1=${IRFS}/tel1/ze${ZENITH}" --irf "tel2=${IRFS}/tel2/ze${ZENITH}" \
  --irf "stereo=${IRFS}/stereo/ze${ZENITH}" --irf "strict_stereo=${IRFS}/strict_stereo/ze${ZENITH}" \
  --dl2 "tel1=${DL2}/tel1" --dl2 "tel2=${DL2}/tel2" --dl2 "stereo=${DL2}/stereo" \
  -o "${IRFS}/comparison/ze${ZENITH}" --title "zenith ${ZENITH} deg, ${IRF_CUTS} cuts"
