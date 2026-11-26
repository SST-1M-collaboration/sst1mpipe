#!/bin/bash

#######################################
#  Merge the DL1 files of simulations #
#######################################
# Job array, one task per merged file: 5 per zenith angle of MC_ZENITHS (task = 5 * zenith index + file).
# For each zenith angle, the first MC_TRAIN_FRACTION of the diffuse gamma and proton files are for the
# training, the others for the tests (IRFs):
#   merged/ze<zenith>/gamma_point.dl1.h5, {gamma_diffuse,proton}_{train,test}.dl1.h5
# The training files of several zenith angles are merged after by mc_merge_train.bash.
# The missing DL1 files (failed jobs) are skipped.

#SBATCH --job-name=mc-merge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=8GB
#SBATCH --output=/home/%u/logs/%x-%A_%a.out
#SBATCH --error=/home/%u/logs/%x-%A_%a.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"

read -ra ZENITHS <<< "${MC_ZENITHS}"
TASK=${SLURM_ARRAY_TASK_ID:-0}
ZENITH=${ZENITHS[TASK / 5]}
NAME=(gamma_point gamma_diffuse_train gamma_diffuse_test proton_train proton_test)
NAME=${NAME[TASK % 5]}
PARTICLE=${NAME%_train}
PARTICLE=${PARTICLE%_test}

# DL1 files of the simtel files of the list of the particle at the zenith angle
mapfile -t FILES < <(sed -e "s|.*/||; s|\.simtel\.gz$|.dl1.h5|; s|^|${MC_DIRECTORY}/dl1/${PARTICLE}/ze${ZENITH}/|" \
  "${MC_DIRECTORY}/lists/${PARTICLE}_ze${ZENITH}.txt")
# training: the first MC_TRAIN_FRACTION of the files, tests: the others
N_TRAIN=$(awk -v n="${#FILES[@]}" -v f="${MC_TRAIN_FRACTION}" 'BEGIN {print int(n * f + 0.5)}')
case "${NAME}" in
  *_train) FILES=("${FILES[@]:0:N_TRAIN}") ;;
  *_test) FILES=("${FILES[@]:N_TRAIN}") ;;
esac

EXISTING=()
for FILE in "${FILES[@]}"; do
  if [[ -f "${FILE}" ]]; then EXISTING+=("${FILE}"); else echo "Missing: ${FILE}" >&2; fi
done
OUTPUT="${MC_DIRECTORY}/merged/ze${ZENITH}/${NAME}.dl1.h5"
mkdir -p "$(dirname "${OUTPUT}")"
echo "${OUTPUT}: ${#EXISTING[@]} files"
ctapipe-merge "${EXISTING[@]}" --output "${OUTPUT}" --overwrite --MergeTool.skip_broken_files=True
