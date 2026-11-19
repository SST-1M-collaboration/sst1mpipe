#!/bin/bash

#####################################################
#  Simulations R0 -> IRFs: submits the jobs (camk)  #
#####################################################
# R0 -> DL1 (one job array per particle and zenith angle) -> merge -> training (all the zenith angles)
# -> DL1 -> DL2 -> IRFs (per zenith angle)
#   bash scripts/camk/submit_mc.bash
# Paths and settings: scripts/camk/sst1mpipe.env (or SST1MPIPE_ENV=<file>).
# MC_START: first step, dl1 (default), merge, train, dl2 or irf: the steps before were done by a
# previous submission (e.g. MC_START=irf IRF_CUTS=paper: the IRFs with other cuts).

export SST1MPIPE_ROOT="${SST1MPIPE_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
CONDA_ENV="" source "${SST1MPIPE_ROOT}/scripts/camk/common.bash"
SCRIPTS="${SST1MPIPE_ROOT}/scripts/camk"
mkdir -p "${MC_DIRECTORY}/lists"

STEPS=(dl1 merge train dl2 irf)
MC_START="${MC_START:-dl1}"
[[ " ${STEPS[*]} " == *" ${MC_START} "* ]] || { echo "MC_START: one of ${STEPS[*]}" >&2; exit 1; }
# run STEP: whether the step is MC_START or after
run() {
  local step
  for step in "${STEPS[@]}"; do
    [[ "${step}" == "${MC_START}" ]] && return 0
    [[ "${step}" == "$1" ]] && return 1
  done
}
# sbatch with the dependency on the previous job (if any)
submit() {
  local dependency=$1
  shift
  if [[ -n "${dependency}" ]]; then
    sbatch --parsable --dependency="${dependency}" "$@"
  else
    sbatch --parsable "$@"
  fi
}

R0_JOBS=()
if run dl1; then
  for ZENITH in ${MC_ZENITHS}; do
    for PARTICLE in gamma_diffuse gamma_point proton; do
      LIST="${MC_DIRECTORY}/lists/${PARTICLE}_ze${ZENITH}.txt"
      PATTERN=$(mc_pattern "${PARTICLE}" "${ZENITH}")
      find "$(dirname "${PATTERN}")" -maxdepth 1 -name "$(basename "${PATTERN}")" 2>/dev/null | sort > "${LIST}" || true
      N=$(wc -l < "${LIST}")
      if (( N == 0 )); then
        echo "No ${PARTICLE} file at zenith ${ZENITH} deg: $(mc_pattern "${PARTICLE}" "${ZENITH}")" >&2
        exit 1
      fi
      N_TASKS=$(( (N + MC_FILES_PER_TASK - 1) / MC_FILES_PER_TASK ))
      R0_JOBS+=("$(sbatch --parsable --array="0-$((N_TASKS - 1))" --export="ALL,PARTICLE=${PARTICLE},ZENITH=${ZENITH}" \
        "${SCRIPTS}/mc_r0_dl1.bash")")
      echo "${PARTICLE} zenith ${ZENITH} deg: ${N} files, ${N_TASKS} tasks, job ${R0_JOBS[-1]}"
    done
  done
fi

read -ra ZENITHS <<< "${MC_ZENITHS}"
N_ZENITHS=${#ZENITHS[@]}
DEPENDENCY=""
# failed DL1 files are skipped by the merge (afterany)
(( ${#R0_JOBS[@]} > 0 )) && DEPENDENCY="afterany:$(join_ids "${R0_JOBS[@]}")"
if run merge; then
  MERGE=$(submit "${DEPENDENCY}" "${SCRIPTS}/mc_merge.bash")
  DEPENDENCY="afterok:${MERGE}"
  echo "merge: job ${MERGE}"
fi
if run train; then
  TRAIN=$(submit "${DEPENDENCY}" --array=0-2 "${SCRIPTS}/mc_train.bash")
  DEPENDENCY="afterok:${TRAIN}"
  echo "training: job ${TRAIN}"
fi
if run dl2; then
  # 3 versions (stereo, tel1, tel2) x 3 files per zenith angle
  DL2=$(submit "${DEPENDENCY}" --array="0-$((9 * N_ZENITHS - 1))" "${SCRIPTS}/mc_dl1_dl2.bash")
  DEPENDENCY="afterok:${DL2}"
  echo "DL2: job ${DL2}"
fi
# 4 modes (tel1, tel2, stereo, strict_stereo) per zenith angle
IRF=$(submit "${DEPENDENCY}" --array="0-$((4 * N_ZENITHS - 1))" "${SCRIPTS}/mc_irf.bash")
echo "IRFs (${IRF_CUTS} cuts): job ${IRF}, outputs in ${MC_DIRECTORY}"
