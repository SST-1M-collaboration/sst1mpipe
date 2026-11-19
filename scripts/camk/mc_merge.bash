#!/bin/bash

#######################################
#  Merge the DL1 files of simulations #
#######################################
# For each zenith angle, the first MC_TRAIN_FRACTION of the diffuse gamma and proton files are for the
# training, the others for the tests (IRFs):
# - merged/{gamma_diffuse,proton}_train.dl1.h5: the training files of all the zenith angles
# - merged/ze<zenith>/{gamma_diffuse,proton}_test.dl1.h5, merged/ze<zenith>/gamma_point.dl1.h5
# The missing DL1 files (failed jobs) are skipped.

#SBATCH --job-name=mc-merge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-04:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=8GB
#SBATCH --output=/home/%u/logs/%x-%j.out
#SBATCH --error=/home/%u/logs/%x-%j.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"

# merge OUTPUT DL1_FILES...: the existing ones
merge() {
  local output=$1 files=()
  shift
  for file in "$@"; do
    if [[ -f "${file}" ]]; then files+=("${file}"); else echo "Missing: ${file}" >&2; fi
  done
  mkdir -p "$(dirname "${output}")"
  echo "${output}: ${#files[@]} files"
  ctapipe-merge "${files[@]}" --output "${output}" --overwrite --MergeTool.skip_broken_files=True
}

# DL1 files of the simtel files of the list of PARTICLE at ZENITH
dl1_files() {
  sed -e "s|.*/||; s|\.simtel\.gz$|.dl1.h5|; s|^|${MC_DIRECTORY}/dl1/$1/ze$2/|" "${MC_DIRECTORY}/lists/$1_ze$2.txt"
}

declare -A TRAIN=([gamma_diffuse]="" [proton]="")
for ZENITH in ${MC_ZENITHS}; do
  mapfile -t FILES < <(dl1_files gamma_point "${ZENITH}")
  merge "${MC_DIRECTORY}/merged/ze${ZENITH}/gamma_point.dl1.h5" "${FILES[@]}"
  for PARTICLE in gamma_diffuse proton; do
    mapfile -t FILES < <(dl1_files "${PARTICLE}" "${ZENITH}")
    N_TRAIN=$(awk -v n="${#FILES[@]}" -v f="${MC_TRAIN_FRACTION}" 'BEGIN {print int(n * f + 0.5)}')
    TRAIN[${PARTICLE}]+=" ${FILES[*]:0:N_TRAIN}"
    merge "${MC_DIRECTORY}/merged/ze${ZENITH}/${PARTICLE}_test.dl1.h5" "${FILES[@]:N_TRAIN}"
  done
done
for PARTICLE in gamma_diffuse proton; do
  # shellcheck disable=SC2086  # file names without spaces
  merge "${MC_DIRECTORY}/merged/${PARTICLE}_train.dl1.h5" ${TRAIN[${PARTICLE}]}
done
