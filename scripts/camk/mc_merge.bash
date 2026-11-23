#!/bin/bash

#######################################
#  Merge the DL1 files of simulations #
#######################################
# For each zenith angle of MC_ZENITHS, the first MC_TRAIN_FRACTION of the diffuse gamma and proton files
# are for the training, the others for the tests (IRFs):
# - merged/ze<zenith>/{gamma_diffuse,proton}_train.dl1.h5, {gamma_diffuse,proton}_test.dl1.h5, gamma_point.dl1.h5
# - merged/train_<MODELS_NAME>/{gamma_diffuse,proton}_train.dl1.h5: the training files of the zenith
#   angles MC_TRAIN_ZENITHS (several: those of the zenith angles merged before are used)
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

for ZENITH in ${MC_ZENITHS}; do
  mapfile -t FILES < <(dl1_files gamma_point "${ZENITH}")
  merge "${MC_DIRECTORY}/merged/ze${ZENITH}/gamma_point.dl1.h5" "${FILES[@]}"
  for PARTICLE in gamma_diffuse proton; do
    mapfile -t FILES < <(dl1_files "${PARTICLE}" "${ZENITH}")
    N_TRAIN=$(awk -v n="${#FILES[@]}" -v f="${MC_TRAIN_FRACTION}" 'BEGIN {print int(n * f + 0.5)}')
    merge "${MC_DIRECTORY}/merged/ze${ZENITH}/${PARTICLE}_train.dl1.h5" "${FILES[@]:0:N_TRAIN}"
    merge "${MC_DIRECTORY}/merged/ze${ZENITH}/${PARTICLE}_test.dl1.h5" "${FILES[@]:N_TRAIN}"
  done
done
# training files of several zenith angles (one zenith angle: its own, see train_directory). Their
# obs_ids (CORSIKA runs) must differ: ctapipe-merge would skip a file with obs_ids already merged.
read -ra TRAIN_ZENITHS <<< "${MC_TRAIN_ZENITHS}"
if (( ${#TRAIN_ZENITHS[@]} > 1 )); then
  for PARTICLE in gamma_diffuse proton; do
    FILES=()
    for ZENITH in "${TRAIN_ZENITHS[@]}"; do
      FILE="${MC_DIRECTORY}/merged/ze${ZENITH}/${PARTICLE}_train.dl1.h5"
      if [[ ! -f "${FILE}" ]]; then
        echo "Missing ${FILE}: merge the zenith angle ${ZENITH} (MC_ZENITHS=\"${ZENITH}\" MC_START=merge)" >&2
        exit 1
      fi
      FILES+=("${FILE}")
    done
    python - "${FILES[@]}" <<'PYTHON'
import sys
from collections import Counter
from ctapipe.io import read_table
obs_ids = {path: set(read_table(path, "/configuration/observation/observation_block")["obs_id"]) for path in sys.argv[1:]}
duplicated = {obs_id for obs_id, n in Counter(o for ids in obs_ids.values() for o in ids).items() if n > 1}
if duplicated:
    files = [path for path, ids in obs_ids.items() if ids & duplicated]
    sys.exit(f"Same obs_ids (CORSIKA runs) in {files}, e.g. {sorted(duplicated)[:5]}: they cannot be merged")
PYTHON
    mkdir -p "$(train_directory)"
    echo "$(train_directory)/${PARTICLE}_train.dl1.h5: ${#FILES[@]} files"
    ctapipe-merge "${FILES[@]}" --output "$(train_directory)/${PARTICLE}_train.dl1.h5" --overwrite
  done
fi
