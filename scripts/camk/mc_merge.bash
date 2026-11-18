#!/bin/bash

#######################################
#  Merge the DL1 files of simulations #
#######################################
# merged/gamma_point.dl1.h5: all the point gammas (performance, point-like IRFs)
# merged/{gamma_diffuse,proton}_train.dl1.h5: the first MC_TRAIN_FRACTION of the files (training)
# merged/{gamma_diffuse,proton}_test.dl1.h5: the others (IRFs)
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
mkdir -p "${MC_DIRECTORY}/merged"

# merge OUTPUT DL1_FILES...: the existing ones
merge() {
  local output=$1 files=()
  shift
  for file in "$@"; do
    if [[ -f "${file}" ]]; then files+=("${file}"); else echo "Missing: ${file}" >&2; fi
  done
  echo "${output}: ${#files[@]} files"
  ctapipe-merge "${files[@]}" --output "${output}" --overwrite --MergeTool.skip_broken_files=True
}

for PARTICLE in gamma_diffuse gamma_point proton; do
  mapfile -t DL1_FILES < <(sed -e "s|.*/||; s|\.simtel\.gz$|.dl1.h5|; s|^|${MC_DIRECTORY}/dl1/${PARTICLE}/|" \
    "${MC_DIRECTORY}/lists/${PARTICLE}.txt")
  if [[ "${PARTICLE}" == gamma_point ]]; then
    merge "${MC_DIRECTORY}/merged/${PARTICLE}.dl1.h5" "${DL1_FILES[@]}"
  else
    N_TRAIN=$(awk -v n="${#DL1_FILES[@]}" -v f="${MC_TRAIN_FRACTION}" 'BEGIN {print int(n * f + 0.5)}')
    merge "${MC_DIRECTORY}/merged/${PARTICLE}_train.dl1.h5" "${DL1_FILES[@]:0:N_TRAIN}"
    merge "${MC_DIRECTORY}/merged/${PARTICLE}_test.dl1.h5" "${DL1_FILES[@]:N_TRAIN}"
  fi
done
