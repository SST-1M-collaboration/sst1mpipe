#!/bin/bash

#############################
#  DL2 -> DL3 of the data   #
#############################
# Every DL2 file of dl2/<cs1, cs2, stereo>/ with the cuts and IRFs of DATA_IRF_DIRECTORY/<mode> of its
# analysis (DATA_IRF_MODES: cs1 -> tel1, cs2 -> tel2, stereo -> stereo), zenith angle closest to
# each observation (ze<zenith>), N_PARALLEL files at a time (sst1mpipe-dl2-dl3): dl3/<cs1, cs2, stereo>/,
# then the index files of gammapy of each directory.

#SBATCH --job-name=data-dl2-dl3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-04:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=/home/%u/logs/%x-%j.out
#SBATCH --error=/home/%u/logs/%x-%j.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"

# dl2/<analysis>/<name>.h5 -> dl3/<analysis>/SST1M_<obs_id>.dl3.fits.gz, IRFs of the mode $2
dl2_dl3() {
  local output
  output="$(dirname "${1/\/dl2\//\/dl3\/}")"
  mkdir -p "${output}"
  sst1mpipe-dl2-dl3 -i "$1" --irf-directory "${DATA_IRF_DIRECTORY}/$2" -o "${output}" --overwrite --no-index
}
export -f dl2_dl3
export DATA_IRF_DIRECTORY

for MAPPING in ${DATA_IRF_MODES}; do
  ANALYSIS="${MAPPING%%=*}"
  MODE="${MAPPING#*=}"
  [[ -d "${DATA_DIRECTORY}/dl2/${ANALYSIS}" ]] || continue
  if [[ ! -d "${DATA_IRF_DIRECTORY}/${MODE}" ]]; then
    echo "${ANALYSIS}: no IRFs ${DATA_IRF_DIRECTORY}/${MODE}, skipped" >&2
    continue
  fi
  find "${DATA_DIRECTORY}/dl2/${ANALYSIS}" -name "*.h5" | sort \
    | xargs -r -P "${N_PARALLEL}" -I{} bash -c 'dl2_dl3 "$1" "$2"' _ {} "${MODE}"
  DL3="${DATA_DIRECTORY}/dl3/${ANALYSIS}"
  if compgen -G "${DL3}/*.dl3.fits.gz" > /dev/null; then
    create_hdu_indexes --input-files "${DL3}/*.dl3.fits.gz" --output-dir "${DL3}"
  fi
done
