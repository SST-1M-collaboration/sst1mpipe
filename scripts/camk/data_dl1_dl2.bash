#!/bin/bash

#############################
#  DL1 -> DL2 of the data   #
#############################
# Every DL1 file (cs1, cs2, stereo) with the models of the node of MODELS_DIRECTORY closest to its
# pointing (sst1mpipe-apply-models), N_PARALLEL files at a time: dl2/<cs1, cs2, stereo>/.

#SBATCH --job-name=data-dl1-dl2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-08:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=/home/%u/logs/%x-%j.out
#SBATCH --error=/home/%u/logs/%x-%j.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"

# dl1/<subdirectory>/<name>.h5 -> dl2/<subdirectory>/<name with dl1 -> dl2>.h5
dl1_dl2() {
  local output="${1/\/dl1\//\/dl2\/}"
  output="$(dirname "${output}")/$(basename "${output//dl1/dl2}")"
  mkdir -p "$(dirname "${output}")"
  sst1mpipe-apply-models -i "$1" -m "${MODELS_DIRECTORY}" -o "${output}" --overwrite --no-progress
}
export -f dl1_dl2
export MODELS_DIRECTORY

find "${DATA_DIRECTORY}/dl1" -name "*.h5" | sort | xargs -r -P "${N_PARALLEL}" -I{} bash -c 'dl1_dl2 "$1"' _ {}
