#!/bin/bash

#####################################################
#  Simulations R0 -> IRFs: submits the jobs (camk)  #
#####################################################
# R0 -> DL1 (one job array per particle) -> merge -> training -> DL1 -> DL2 -> IRFs
#   bash scripts/camk/submit_mc.bash
# Paths and settings: scripts/camk/sst1mpipe.env (or SST1MPIPE_ENV=<file>).

export SST1MPIPE_ROOT="${SST1MPIPE_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
CONDA_ENV="" source "${SST1MPIPE_ROOT}/scripts/camk/common.bash"
SCRIPTS="${SST1MPIPE_ROOT}/scripts/camk"
mkdir -p "${MC_DIRECTORY}/lists"

declare -A INPUTS=([gamma_diffuse]="${MC_GAMMA_DIFFUSE}" [gamma_point]="${MC_GAMMA_POINT}" [proton]="${MC_PROTON}")
R0_JOBS=()
for PARTICLE in gamma_diffuse gamma_point proton; do
  LIST="${MC_DIRECTORY}/lists/${PARTICLE}.txt"
  # shellcheck disable=SC2086  # glob pattern
  ls -1 ${INPUTS[${PARTICLE}]} > "${LIST}"
  N=$(wc -l < "${LIST}")
  R0_JOBS+=("$(sbatch --parsable --array="0-$((N - 1))" --export="ALL,PARTICLE=${PARTICLE}" "${SCRIPTS}/mc_r0_dl1.bash")")
  echo "${PARTICLE}: ${N} files, job ${R0_JOBS[-1]}"
done

# failed DL1 files are skipped by the merge (afterany)
MERGE=$(sbatch --parsable --dependency="afterany:$(join_ids "${R0_JOBS[@]}")" "${SCRIPTS}/mc_merge.bash")
TRAIN=$(sbatch --parsable --array=0-2 --dependency="afterok:${MERGE}" "${SCRIPTS}/mc_train.bash")
DL2=$(sbatch --parsable --array=0-2 --dependency="afterok:${TRAIN}" "${SCRIPTS}/mc_dl1_dl2.bash")
IRF=$(sbatch --parsable --dependency="afterok:${DL2}" "${SCRIPTS}/mc_irf.bash")
echo "merge ${MERGE}, training ${TRAIN}, DL2 ${DL2}, IRFs ${IRF}: outputs in ${MC_DIRECTORY}"
