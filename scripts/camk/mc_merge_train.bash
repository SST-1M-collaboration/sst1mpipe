#!/bin/bash

###############################################################
#  Merge the training files of several zenith angles          #
###############################################################
# Job array, one task per particle (0: gamma_diffuse, 1: proton): the training files
# merged/ze<zenith>/<particle>_train.dl1.h5 of the zenith angles of MC_TRAIN_ZENITHS (merged by
# mc_merge.bash, in this submission or before) in merged/train_<MODELS_NAME>/<particle>_train.dl1.h5.
# Their obs_ids (CORSIKA runs) must differ: ctapipe-merge would skip a file with obs_ids already merged.

#SBATCH --job-name=mc-merge-train
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=8GB
#SBATCH --output=/home/%u/logs/%x-%A_%a.out
#SBATCH --error=/home/%u/logs/%x-%A_%a.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"

PARTICLES=(gamma_diffuse proton)
PARTICLE=${PARTICLES[${SLURM_ARRAY_TASK_ID:-0}]}
read -ra TRAIN_ZENITHS <<< "${MC_TRAIN_ZENITHS}"
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
OUTPUT="$(train_directory)/${PARTICLE}_train.dl1.h5"
mkdir -p "$(dirname "${OUTPUT}")"
echo "${OUTPUT}: ${#FILES[@]} files"
ctapipe-merge "${FILES[@]}" --output "${OUTPUT}" --overwrite
