# Sourced by the scripts of scripts/camk: strict mode, environment variables, conda environment.
# SST1MPIPE_ROOT: sst1mpipe directory, exported by the submit scripts (sbatch runs a copy of the
# job scripts: they cannot find their directory). SST1MPIPE_ENV: another environment file.

set -euo pipefail
export SST1MPIPE_ROOT="${SST1MPIPE_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
source "${SST1MPIPE_ENV:-${SST1MPIPE_ROOT}/scripts/camk/sst1mpipe.env}"
mkdir -p ~/logs
# Python traceback in the log of the job when a compiled library crashes (segmentation fault)
export PYTHONFAULTHANDLER=1

if [[ -n "${CONDA_ENV}" ]]; then
  # a job is a new shell: mamba is initialised before "mamba activate", without "set -u"
  set +u
  eval "$(mamba shell hook --shell bash)"
  mamba activate "${CONDA_ENV}"
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

# -r <model> of the models (*.pkl) of a directory, for sst1mpipe-apply-models
models_args() {
  local model
  for model in "$1"/*.pkl; do
    printf -- '-r\n%s\n' "${model}"
  done
}

# simtel files of a particle (MC_GAMMA_DIFFUSE, MC_GAMMA_POINT, MC_PROTON) at a zenith angle
mc_pattern() {
  local -A patterns=([gamma_diffuse]="${MC_GAMMA_DIFFUSE}" [gamma_point]="${MC_GAMMA_POINT}" [proton]="${MC_PROTON}")
  local pattern="${patterns[$1]}"
  echo "${pattern//%zenith%/$2}"
}
