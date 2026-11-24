#!/bin/bash

#############################
#  DL1 -> DL2 of the data   #
#############################
# Every DL1 file with the models of MODELS_DIRECTORY (sst1mpipe-apply-models): cs1, cs2 with the mono
# models (MODELS_DIRECTORY/mono), stereo with the stereo models (MODELS_DIRECTORY/stereo),
# N_PARALLEL files at a time, DL2_N_JOBS threads each: dl2/<cs1, cs2, stereo>/.

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

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}/scripts/daily/common.bash"

# dl1/<subdirectory>/<name>.h5 -> dl2/<subdirectory>/<name with dl1 -> dl2>.h5
dl1_dl2() {
  local output="${1/\/dl1\//\/dl2\/}"
  output="$(dirname "${output}")/$(basename "${output//dl1/dl2}")"
  mkdir -p "$(dirname "${output}")"
  local models set=mono
  [[ "$1" == */dl1/stereo/* ]] && set=stereo
  mapfile -t models < <(models_args "${MODELS_DIRECTORY}/${set}")
  (( ${#models[@]} > 0 )) || return 1
  sst1mpipe-apply-models -i "$1" "${models[@]}" -o "${output}" --n-jobs "${DL2_N_JOBS}" --overwrite --no-progress
}
export -f dl1_dl2 models_args
export MODELS_DIRECTORY DL2_N_JOBS

# not the night files merged for the stereo (dl1/merged)
find "${DATA_DIRECTORY}/dl1" -path "*/dl1/merged" -prune -o -name "*.h5" -print | sort \
  | xargs -r -P "${N_PARALLEL}" -I{} bash -c 'dl1_dl2 "$1"' _ {}
