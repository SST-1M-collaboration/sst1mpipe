#!/bin/bash

####################################################################
#  Data of a night R0 -> DL3 with SLURM: submits the jobs (camk)   #
####################################################################
# R0 -> DL1 (one job array per telescope) -> DL1 stereo -> DL2 -> DL3, logs in ~/logs
#   DATE=20260120 bash scripts/daily/submit_data.bash
# DATA_STOP: last step, dl1 (monitoring and R0 -> DL1 only), stereo, dl2 or dl3 (default).
# RAW_LIST: only the raw files of this list (see raw_files in common.bash), e.g. the runs of a source.
# Several nights: submit_nights.bash; a dataset (ECSV): scripts/dataset/submit_dataset.bash.
# LAST_JOBS_FILE: file to which the id of the last job of the night is appended (dependency of a later job).
# Settings: scripts/daily/camk.env (DAILY_ENV=<file>: other settings). Without SLURM: run_daily.bash.

export SST1MPIPE_ROOT="${SST1MPIPE_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
export DAILY_ENV="${DAILY_ENV:-${SST1MPIPE_ROOT}/scripts/daily/camk.env}"
CONDA_ENV="" source "${SST1MPIPE_ROOT}/scripts/daily/common.bash"
export DATE
DATA_STOP="${DATA_STOP:-dl3}"
[[ " dl1 stereo dl2 dl3 " == *" ${DATA_STOP} "* ]] || { echo "DATA_STOP: dl1, stereo, dl2 or dl3" >&2; exit 1; }
[[ -z "${RAW_LIST:-}" || -f "${RAW_LIST}" ]] || { echo "RAW_LIST: no file ${RAW_LIST}" >&2; exit 1; }
SCRIPTS="${SST1MPIPE_ROOT}/scripts/daily"
# ids of the last jobs of the night, appended to LAST_JOBS_FILE
record_jobs() {
  [[ -z "${LAST_JOBS_FILE:-}" ]] || printf '%s\n' "$@" >> "${LAST_JOBS_FILE}"
}
mkdir -p "${DATA_DIRECTORY}/lists" ~/logs

R0_JOBS=()
for TEL in 1 2; do
  LIST="${DATA_DIRECTORY}/lists/cs${TEL}.txt"
  raw_files "${TEL}" > "${LIST}"
  N=$(wc -l < "${LIST}")
  if (( N == 0 )); then
    echo "cs${TEL}: no raw file of ${DATE}" >&2
    continue
  fi
  R0_JOBS+=("$(sbatch --parsable --array="0-$((N - 1))" --export="ALL,TEL=${TEL}" "${SCRIPTS}/data_r0_dl1.bash")")
  echo "cs${TEL}: ${N} files, job ${R0_JOBS[-1]}"
done
(( ${#R0_JOBS[@]} > 0 )) || exit 1
if [[ "${DATA_STOP}" == dl1 ]]; then
  record_jobs "${R0_JOBS[@]}"
  echo "outputs in ${DATA_DIRECTORY}/{monitoring,dl1}"
  exit 0
fi

# failed DL1 files are skipped (afterany)
DEPENDENCY="afterany:$(join_ids "${R0_JOBS[@]}")"
if (( ${#R0_JOBS[@]} == 2 )); then
  STEREO=$(sbatch --parsable --dependency="${DEPENDENCY}" "${SCRIPTS}/data_dl1_stereo.bash")
  DEPENDENCY="afterany:${STEREO}"
  echo "stereo: job ${STEREO}"
fi
if [[ "${DATA_STOP}" == stereo ]]; then
  IFS=: read -ra LAST <<< "${DEPENDENCY#afterany:}"
  record_jobs "${LAST[@]}"
  exit 0
fi
DL2=$(sbatch --parsable --dependency="${DEPENDENCY}" "${SCRIPTS}/data_dl1_dl2.bash")
echo "DL2: job ${DL2}"
if [[ "${DATA_STOP}" == dl2 ]]; then
  record_jobs "${DL2}"
  exit 0
fi
DL3=$(sbatch --parsable --dependency="afterany:${DL2}" "${SCRIPTS}/data_dl2_dl3.bash")
echo "DL3: job ${DL3}, outputs in ${DATA_DIRECTORY}"
record_jobs "${DL3}"
