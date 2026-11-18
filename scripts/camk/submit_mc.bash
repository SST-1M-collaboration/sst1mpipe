#!/bin/bash

#####################################################
#  Simulations R0 -> IRFs: submits the jobs (camk)  #
#####################################################
# R0 -> DL1 (one job array per particle and zenith angle) -> merge -> training (all the zenith angles)
# -> DL1 -> DL2 -> IRFs (per zenith angle)
#   bash scripts/camk/submit_mc.bash
# Paths and settings: scripts/camk/sst1mpipe.env (or SST1MPIPE_ENV=<file>).

export SST1MPIPE_ROOT="${SST1MPIPE_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
CONDA_ENV="" source "${SST1MPIPE_ROOT}/scripts/camk/common.bash"
SCRIPTS="${SST1MPIPE_ROOT}/scripts/camk"
mkdir -p "${MC_DIRECTORY}/lists"

R0_JOBS=()
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

# failed DL1 files are skipped by the merge (afterany)
MERGE=$(sbatch --parsable --dependency="afterany:$(join_ids "${R0_JOBS[@]}")" "${SCRIPTS}/mc_merge.bash")
TRAIN=$(sbatch --parsable --array=0-2 --dependency="afterok:${MERGE}" "${SCRIPTS}/mc_train.bash")
read -ra ZENITHS <<< "${MC_ZENITHS}"
N_ZENITHS=${#ZENITHS[@]}
DL2=$(sbatch --parsable --array="0-$((3 * N_ZENITHS - 1))" --dependency="afterok:${TRAIN}" "${SCRIPTS}/mc_dl1_dl2.bash")
IRF=$(sbatch --parsable --array="0-$((N_ZENITHS - 1))" --dependency="afterok:${DL2}" "${SCRIPTS}/mc_irf.bash")
echo "merge ${MERGE}, training ${TRAIN}, DL2 ${DL2}, IRFs ${IRF}: outputs in ${MC_DIRECTORY}"
