#!/bin/bash

#####################################################################
#  Data of a dataset (ECSV) with SLURM: submits the jobs (camk)      #
#####################################################################
# The runs of a dataset (sst1mpipe.io.dataset, e.g. sst1mpipe/resources/datasets/crab_paper_2025.ecsv),
# night by night with the daily analysis (scripts/daily/submit_data.bash: monitoring and R0 -> DL1 per
# raw file -> DL1 stereo -> DL2 -> DL3) and the models and IRFs of the simulations
# (scripts/camk/submit_mc.bash: MC_DIRECTORY, MODELS_NAME, IRF_CUTS), then the index files of gammapy of
# all the nights (dataset_dl3_index.bash):
#   DATASET=crab_paper_2025 bash scripts/dataset/submit_dataset.bash
#   DATASET=crab_paper_2025 DATA_STOP=dl1 bash scripts/dataset/submit_dataset.bash
#
# DATASET: name of a dataset of sst1mpipe, or an ECSV file. DATASET_ROOT: root of the raw files (default:
# the root of the dataset if it exists, else RAW_DIRECTORY of camk.env). DATASET_NIGHTS: only these nights
# (YYYYMMDD, separated by spaces). DATA_STOP: last step (dl1, stereo, dl2, dl3), see submit_data.bash.
# Outputs: DATASETS_DIRECTORY/<name>/ (default: next to DAILY_DIRECTORY, .../sst1mpipe/datasets/<name>):
#   dataset/   the dataset, the raw files, the configurations, settings and version (reproducibility)
#   <night>/   outputs of the night (lists, monitoring, dl1, dl2, dl3), as the daily analysis
#   dl3/<cs1, cs2, stereo>/   index files of gammapy of the DL3 files of all the nights

export SST1MPIPE_ROOT="${SST1MPIPE_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
export DAILY_ENV="${DAILY_ENV:-${SST1MPIPE_ROOT}/scripts/daily/camk.env}"
CONDA_ENV="" source "${SST1MPIPE_ROOT}/scripts/daily/common.bash"
DATA_STOP="${DATA_STOP:-dl3}"
[[ " dl1 stereo dl2 dl3 " == *" ${DATA_STOP} "* ]] || { echo "DATA_STOP: dl1, stereo, dl2 or dl3" >&2; exit 1; }

# the dataset: a file, or a dataset of sst1mpipe by name
DATASET="${DATASET:?DATASET: name of a dataset of sst1mpipe (sst1mpipe/resources/datasets) or ECSV file}"
if [[ -f "${DATASET}" ]]; then
  DATASET_FILE="$(realpath "${DATASET}")"
else
  DATASET_FILE="${SST1MPIPE_ROOT}/sst1mpipe/resources/datasets/${DATASET%.ecsv}.ecsv"
  [[ -f "${DATASET_FILE}" ]] || { echo "No dataset ${DATASET} (file, or ${DATASET_FILE})" >&2; exit 1; }
fi
if head -c 40 "${DATASET_FILE}" | grep -q '^version https://git-lfs'; then
  echo "${DATASET_FILE} is a git LFS pointer: git lfs pull -I 'sst1mpipe/resources/datasets/*'" >&2
  exit 1
fi
# header (sst1mpipe.io.dataset): name, root of the raw files
meta() {
  sed -n '/^# meta:/,/^# schema:/p' "${DATASET_FILE}" | sed -n "s/^# - {$1: \(.*\)}$/\1/p" | head -1 \
    | sed "s/^'\(.*\)'$/\1/"
}
NAME="$(meta name)"
NAME="${NAME:-$(basename "${DATASET_FILE}" .ecsv)}"
if [[ -z "${DATASET_ROOT:-}" ]]; then
  DATASET_ROOT="$(meta root)"
  [[ -d "${DATASET_ROOT}" ]] || DATASET_ROOT="${RAW_DIRECTORY}"
fi
DATASETS_DIRECTORY="${DATASETS_DIRECTORY:-$(dirname "${DAILY_DIRECTORY}")/datasets}"
OUTPUT="${DATASETS_DIRECTORY}/${NAME}"
mkdir -p "${OUTPUT}/dataset" ~/logs

# raw files of the dataset (absolute paths): the RAW_LIST of the nights
RAW_FILES="${OUTPUT}/dataset/raw_files.txt"
grep -v '^#' "${DATASET_FILE}" | tail -n +2 | awk -v root="${DATASET_ROOT}" 'NF >= 3 {print root "/" $3}' > "${RAW_FILES}"
N_FILES=$(wc -l < "${RAW_FILES}")
N_MISSING=0
while IFS= read -r FILE; do
  [[ -f "${FILE}" ]] || N_MISSING=$((N_MISSING + 1))
done < "${RAW_FILES}"
echo "dataset ${NAME}: ${N_FILES} raw files in ${DATASET_ROOT} (${N_MISSING} missing)"
(( N_MISSING < N_FILES )) || { echo "No raw file found: DATASET_ROOT?" >&2; exit 1; }

if [[ -n "${DATASET_NIGHTS:-}" ]]; then
  read -ra NIGHT_LIST <<< "${DATASET_NIGHTS}"
else
  mapfile -t NIGHT_LIST < <(list_nights "${RAW_FILES}")
fi

# models and IRFs of the simulations (scripts/camk/submit_mc.bash)
if [[ " dl2 dl3 " == *" ${DATA_STOP} "* ]]; then
  for SET in mono stereo; do
    compgen -G "${MODELS_DIRECTORY}/${SET}/*.pkl" > /dev/null \
      || { echo "No models in ${MODELS_DIRECTORY}/${SET}: scripts/camk/submit_mc.bash (MODELS_NAME)" >&2; exit 1; }
  done
fi
if [[ "${DATA_STOP}" == dl3 ]]; then
  for MAPPING in ${DATA_IRF_MODES}; do
    [[ -d "${DATA_IRF_DIRECTORY}/${MAPPING#*=}" ]] \
      || echo "warning: no IRFs ${DATA_IRF_DIRECTORY}/${MAPPING#*=} (${MAPPING%%=*} without DL3)" >&2
  done
fi

# reproducibility: the dataset, the configurations, the settings and the version
cp "${DATASET_FILE}" "${OUTPUT}/dataset/"
cp "${DATA_CONFIG}" "${MONITORING_CONFIG}" "${OUTPUT}/dataset/"
{
  echo "submitted: $(date -u +%Y-%m-%dT%H:%M:%SZ) by ${USER:-$(id -un)} on $(hostname)"
  echo "sst1mpipe: $(git -C "${SST1MPIPE_ROOT}" describe --always --dirty 2>/dev/null || echo unknown)" \
    "($(git -C "${SST1MPIPE_ROOT}" rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown))"
  for VARIABLE in DATASET_FILE DATASET_ROOT DATA_STOP DATA_CONFIG MONITORING_CONFIG MC_DIRECTORY MODELS_NAME \
      IRF_CUTS MODELS_DIRECTORY DATA_IRF_DIRECTORY DATA_IRF_MODES; do
    echo "${VARIABLE}=${!VARIABLE}"
  done
  echo "nights: ${NIGHT_LIST[*]}"
} > "${OUTPUT}/dataset/settings.txt"

# the nights, as the daily analysis
LAST_JOBS="${OUTPUT}/dataset/last_jobs.txt"
: > "${LAST_JOBS}"
FAILED=()
for NIGHT in "${NIGHT_LIST[@]}"; do
  echo "=== ${NIGHT}"
  env -u DATA_DIRECTORY DATE="${NIGHT}" DAILY_DIRECTORY="${OUTPUT}" RAW_LIST="${RAW_FILES}" DATA_STOP="${DATA_STOP}" \
    LAST_JOBS_FILE="${LAST_JOBS}" bash "${SST1MPIPE_ROOT}/scripts/daily/submit_data.bash" || FAILED+=("${NIGHT}")
done
(( ${#FAILED[@]} == 0 )) || echo "nights not submitted (no raw file or error): ${FAILED[*]}" >&2

# index files of all the nights
if [[ "${DATA_STOP}" == dl3 && -s "${LAST_JOBS}" ]]; then
  mapfile -t JOBS < "${LAST_JOBS}"
  INDEX=$(sbatch --parsable --dependency="afterany:$(join_ids "${JOBS[@]}")" \
    --export="ALL,DATASET_DIRECTORY=${OUTPUT}" "${SST1MPIPE_ROOT}/scripts/dataset/dataset_dl3_index.bash")
  echo "index files of the dataset: job ${INDEX}, ${OUTPUT}/dl3/<cs1, cs2, stereo>"
fi
echo "outputs in ${OUTPUT}"
