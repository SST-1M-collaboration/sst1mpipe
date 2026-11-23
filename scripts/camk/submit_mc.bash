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

read -ra ZENITHS <<< "${MC_ZENITHS}"
N_ZENITHS=${#ZENITHS[@]}

R0_JOBS=()
if run dl1; then
  # lists of the simtel files of each particle and zenith angle
  LISTS=()
  for ZENITH in ${MC_ZENITHS}; do
    for PARTICLE in gamma_diffuse gamma_point proton; do
      LIST="${MC_DIRECTORY}/lists/${PARTICLE}_ze${ZENITH}.txt"
      PATTERN=$(mc_pattern "${PARTICLE}" "${ZENITH}")
      find "$(dirname "${PATTERN}")" -maxdepth 1 -name "$(basename "${PATTERN}")" 2>/dev/null | sort > "${LIST}" || true
      if [[ ! -s "${LIST}" ]]; then
        echo "No ${PARTICLE} file at zenith ${ZENITH} deg: ${PATTERN}" >&2
        exit 1
      fi
      LISTS+=("${PARTICLE}:${ZENITH}:${LIST}")
    done
  done
  # check, before submitting anything: the array tasks and the jobs of the next steps (merge, training,
  # DL2, IRFs, comparison) with the jobs already in the queue within MC_MAX_SUBMIT
  n_tasks() {  # array tasks of all the lists with $1 files per task
    local total=0 entry n
    for entry in "${LISTS[@]}"; do
      n=$(wc -l < "${entry##*:}")
      total=$(( total + (n + $1 - 1) / $1 ))
    done
    echo "${total}"
  }
  QUEUED=$(command -v squeue > /dev/null && squeue -u "${USER:-$(id -un)}" -h -r 2>/dev/null | wc -l || echo 0)
  AVAILABLE=$(( MC_MAX_SUBMIT - QUEUED - 1 - 3 - 14 * N_ZENITHS ))
  TOTAL=$(n_tasks "${MC_FILES_PER_TASK}")
  if (( TOTAL > AVAILABLE )); then
    MINIMUM=${MC_FILES_PER_TASK}
    while (( AVAILABLE > 0 && $(n_tasks "${MINIMUM}") > AVAILABLE )); do MINIMUM=$(( MINIMUM + 1 )); done
    echo "MC_FILES_PER_TASK=${MC_FILES_PER_TASK}: ${TOTAL} array tasks, ${AVAILABLE} jobs available (MC_MAX_SUBMIT" \
      "${MC_MAX_SUBMIT}, ${QUEUED} jobs in the queue, $(( 4 + 14 * N_ZENITHS )) for the next steps)." \
      "Nothing submitted: MC_FILES_PER_TASK=${MINIMUM} or more fits" >&2
    exit 1
  fi
  echo "${MC_FILES_PER_TASK} files per task: ${TOTAL} array tasks (${QUEUED} jobs in the queue, limit ${MC_MAX_SUBMIT})"
  for ENTRY in "${LISTS[@]}"; do
    IFS=: read -r PARTICLE ZENITH LIST <<< "${ENTRY}"
    N=$(wc -l < "${LIST}")
    N_TASKS=$(( (N + MC_FILES_PER_TASK - 1) / MC_FILES_PER_TASK ))
    R0_JOBS+=("$(sbatch --parsable --array="0-$((N_TASKS - 1))" --export="ALL,PARTICLE=${PARTICLE},ZENITH=${ZENITH}" \
      "${SCRIPTS}/mc_r0_dl1.bash")")
    echo "${PARTICLE} zenith ${ZENITH} deg: ${N} files, ${N_TASKS} tasks, job ${R0_JOBS[-1]}"
  done
fi

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
  # starting from the DL2: the models must be there
  if [[ "${MC_START}" == dl2 ]] && ! compgen -G "${MC_DIRECTORY}/models/${MODELS_NAME}/*.pkl" > /dev/null; then
    echo "No models in ${MC_DIRECTORY}/models/${MODELS_NAME} (MODELS_NAME from MC_TRAIN_ZENITHS=\"${MC_TRAIN_ZENITHS}\")" >&2
    exit 1
  fi
  # 3 versions (stereo, tel1, tel2) x 3 files per zenith angle
  DL2=$(submit "${DEPENDENCY}" --array="0-$((9 * N_ZENITHS - 1))" --cpus-per-task="${MC_DL2_CPUS}" \
    --mem-per-cpu="${MC_DL2_MEM_PER_CPU}" "${SCRIPTS}/mc_dl1_dl2.bash")
  DEPENDENCY="afterok:${DL2}"
  echo "DL2: job ${DL2}"
fi
# 4 modes (tel1, tel2, stereo, strict_stereo) per zenith angle
IRF=$(submit "${DEPENDENCY}" --array="0-$((4 * N_ZENITHS - 1))" "${SCRIPTS}/mc_irf.bash")
# figures with the four modes (afterany: the modes that succeeded)
COMPARISON=$(submit "afterany:${IRF}" --array="0-$((N_ZENITHS - 1))" "${SCRIPTS}/mc_irf_comparison.bash")
echo "IRFs (${IRF_CUTS} cuts): job ${IRF}, comparison ${COMPARISON}, outputs in ${MC_DIRECTORY}"
