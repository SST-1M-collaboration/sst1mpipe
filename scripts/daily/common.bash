# Sourced by the scripts of scripts/daily: strict mode, settings, conda environment.
# SST1MPIPE_ROOT: sst1mpipe directory, exported by submit_data.bash / run_daily.bash (sbatch runs a copy of
# the job scripts: they cannot find their directory). DAILY_ENV: file of the settings (default daily.env).

set -euo pipefail
export SST1MPIPE_ROOT="${SST1MPIPE_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
export DAILY_ENV="${DAILY_ENV:-${SST1MPIPE_ROOT}/scripts/daily/daily.env}"
source "${DAILY_ENV}"
# Python traceback in the log when a compiled library crashes (segmentation fault)
export PYTHONFAULTHANDLER=1

# conda environment: conda or mamba from the PATH, or from CONDA_ROOT / the usual installations (cron
# jobs have no initialised shell)
if [[ -n "${CONDA_ENV}" && "${CONDA_DEFAULT_ENV:-}" != "${CONDA_ENV}" ]]; then
  set +u
  for CONDA in conda mamba "${CONDA_ROOT:-}/bin/conda" ~/miniforge3/bin/conda ~/mambaforge/bin/conda \
      ~/miniconda3/bin/conda ~/anaconda3/bin/conda; do
    command -v "${CONDA}" > /dev/null 2>&1 && break
  done
  if [[ "$(basename "${CONDA}")" == mamba ]]; then
    eval "$("${CONDA}" shell hook --shell bash)"
    mamba activate "${CONDA_ENV}"
  else
    eval "$("${CONDA}" shell.bash hook)"
    conda activate "${CONDA_ENV}"
  fi
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

# -r <model> of the models (*.pkl) of a directory, for sst1mpipe-apply-models; error without models
models_args() {
  local model models=()
  for model in "$1"/*.pkl; do
    [[ -f "${model}" ]] && models+=("${model}")
  done
  if (( ${#models[@]} == 0 )); then
    echo "No models (*.pkl) in $1: MODELS_DIRECTORY?" >&2
    return 1
  fi
  for model in "${models[@]}"; do
    printf -- '-r\n%s\n' "${model}"
  done
}

# raw files of the night DATE of the telescope $1 (1 or 2), sorted
raw_files() {
  local directory="${RAW_DIRECTORY}/cs$1/data/raw/${DATE:0:4}/${DATE:4:2}/${DATE:6:2}"
  [[ -d "${directory}" ]] || return 0
  find "${directory}" -name "SST1M$1_*.fits.fz" | sort
}
