#!/bin/bash

######################################################################
#  Daily analysis: data of a night R0 -> DL3 on one machine          #
######################################################################
# Without SLURM, e.g. on Calculus (32 cores), every morning by cron (scripts/daily/crontab):
#   bash scripts/daily/run_daily.bash                # the night before (yesterday)
#   DATE=20260120 bash scripts/daily/run_daily.bash  # a given night
# R0 -> DL1 (N_PARALLEL raw files at a time, default: the cores of the machine) -> DL1 stereo -> DL2 -> DL3,
# in ${DAILY_DIRECTORY}/<DATE>/, logs of each step and raw file in <DATE>/logs/. A step that fails for
# some files does not stop the next steps (the files that succeeded are processed).
# Settings: scripts/daily/daily.env (DAILY_ENV=<file>: other settings). One run at a time (lock file).

export SST1MPIPE_ROOT="${SST1MPIPE_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
source "${SST1MPIPE_ROOT}/scripts/daily/common.bash"
export DATE DATA_DIRECTORY N_PARALLEL
SCRIPTS="${SST1MPIPE_ROOT}/scripts/daily"
LOGS="${DATA_DIRECTORY}/logs"

# one run at a time: the steps use all the cores
mkdir -p "${DAILY_DIRECTORY}"
exec 9> "${DAILY_DIRECTORY}/.run_daily.lock"
if ! flock -n 9; then
  echo "$(date -Is) ${DATE}: another daily analysis is running, stopped" >&2
  exit 1
fi

mkdir -p "${DATA_DIRECTORY}/lists" "${LOGS}"
echo "$(date -Is) ${DATE}: sst1mpipe $(python -c 'import sst1mpipe; print(sst1mpipe.__version__)'), ${N_PARALLEL} files in parallel, outputs in ${DATA_DIRECTORY}"

# step NAME SCRIPT: runs a step, its log in logs/NAME.log; the failure is reported, not fatal
step() {
  local name=$1
  shift
  echo "$(date -Is) ${name}"
  if ! bash "$@" > "${LOGS}/${name}.log" 2>&1; then
    echo "$(date -Is) ${name}: failed (or some files failed), see ${LOGS}/${name}.log" >&2
  fi
}

N_TELESCOPES=0
for TEL in 1 2; do
  LIST="${DATA_DIRECTORY}/lists/cs${TEL}.txt"
  raw_files "${TEL}" > "${LIST}"
  N=$(wc -l < "${LIST}")
  echo "cs${TEL}: ${N} raw files"
  if (( N > 0 )); then N_TELESCOPES=$(( N_TELESCOPES + 1 )); fi
done
if (( N_TELESCOPES == 0 )); then
  echo "$(date -Is) ${DATE}: no raw file in ${RAW_DIRECTORY}, nothing to do"
  exit 0
fi

# R0 -> DL1: one raw file per process, its log in logs/r0_dl1/<run>.log
mkdir -p "${LOGS}/r0_dl1"
echo "$(date -Is) R0 -> DL1"
cat "${DATA_DIRECTORY}"/lists/cs*.txt | xargs -r -P "${N_PARALLEL}" -I{} bash -c \
  'bash "$1" "$2" > "$3/$(basename "$2" .fits.fz).log" 2>&1 || echo "$2: failed, see $3/$(basename "$2" .fits.fz).log" >&2' \
  _ "${SCRIPTS}/data_r0_dl1.bash" {} "${LOGS}/r0_dl1"

if (( N_TELESCOPES == 2 )); then
  step dl1_stereo "${SCRIPTS}/data_dl1_stereo.bash"
fi
step dl1_dl2 "${SCRIPTS}/data_dl1_dl2.bash"
step dl2_dl3 "${SCRIPTS}/data_dl2_dl3.bash"
echo "$(date -Is) ${DATE}: done, DL3 in ${DATA_DIRECTORY}/dl3"
