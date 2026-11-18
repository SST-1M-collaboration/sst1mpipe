# Sourced by the scripts of scripts/camk: strict mode, environment variables, conda environment.
# SST1MPIPE_ROOT: sst1mpipe directory, exported by the submit scripts (sbatch runs a copy of the
# job scripts: they cannot find their directory). SST1MPIPE_ENV: another environment file.
# CONDA_ENV_VARIABLE: variable with the conda environment to activate (CONDA_ENV by default).

set -euo pipefail
export SST1MPIPE_ROOT="${SST1MPIPE_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
source "${SST1MPIPE_ENV:-${SST1MPIPE_ROOT}/scripts/camk/sst1mpipe.env}"
mkdir -p ~/logs

CONDA_ENV_VARIABLE="${CONDA_ENV_VARIABLE:-CONDA_ENV}"
ENVIRONMENT="${!CONDA_ENV_VARIABLE}"
if [[ -n "${ENVIRONMENT}" ]]; then
  # a job is a new shell: mamba is initialised before "mamba activate", without "set -u"
  set +u
  eval "$(mamba shell hook --shell bash)"
  mamba activate "${ENVIRONMENT}"
  set -u
fi

# line SLURM_ARRAY_TASK_ID (0 without job array) of a file list
task_line() {
  sed -n "$(( ${SLURM_ARRAY_TASK_ID:-0} + 1 ))p" "$1"
}

# colon separated job ids, for --dependency
join_ids() {
  local IFS=:
  echo "$*"
}
