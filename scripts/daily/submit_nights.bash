#!/bin/bash

##############################################################
#  Data of several nights with SLURM: submit_data.bash each  #
##############################################################
# The nights of NIGHTS (YYYYMMDD, separated by spaces), or of the raw files of RAW_LIST (directories
# .../YYYY/MM/DD/). Same settings as submit_data.bash (DATA_STOP, RAW_LIST, camk.env), e.g. monitoring
# and R0 -> DL1 of the runs of a list:
#   RAW_LIST=crab_runs.txt DATA_STOP=dl1 bash scripts/daily/submit_nights.bash
#   NIGHTS="20231012 20231013" DATA_STOP=dl1 bash scripts/daily/submit_nights.bash
# A night without raw files (of the list) is skipped.

set -uo pipefail
export SST1MPIPE_ROOT="${SST1MPIPE_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
SCRIPTS="${SST1MPIPE_ROOT}/scripts/daily"
if [[ -n "${RAW_LIST:-}" ]]; then
  [[ -f "${RAW_LIST}" ]] || { echo "RAW_LIST: no file ${RAW_LIST}" >&2; exit 1; }
  export RAW_LIST="$(realpath "${RAW_LIST}")"
fi
if [[ -n "${NIGHTS:-}" ]]; then
  read -ra LIST <<< "${NIGHTS}"
elif [[ -n "${RAW_LIST:-}" ]]; then
  mapfile -t LIST < <(grep -oE '/[0-9]{4}/[0-9]{2}/[0-9]{2}/' "${RAW_LIST}" | tr -d / | sort -u)
else
  echo "NIGHTS or RAW_LIST needed" >&2
  exit 1
fi
(( ${#LIST[@]} > 0 )) || { echo "No night (NIGHTS, or directories YYYY/MM/DD in RAW_LIST)" >&2; exit 1; }

echo "${#LIST[@]} nights: ${LIST[*]}"
FAILED=()
for NIGHT in "${LIST[@]}"; do
  echo "=== ${NIGHT}"
  DATE="${NIGHT}" bash "${SCRIPTS}/submit_data.bash" || FAILED+=("${NIGHT}")
done
if (( ${#FAILED[@]} > 0 )); then
  echo "nights not submitted (no raw file or error): ${FAILED[*]}" >&2
fi
