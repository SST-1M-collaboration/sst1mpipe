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

# raw files of the night DATE of the telescope $1 (1 or 2), sorted: all the files of the night, or with
# RAW_LIST (file with one raw file per line: absolute, or relative to ${RAW_DIRECTORY}/cs<tel>/, e.g.
# data/raw/2023/10/12/SST1M1_20231012_0123.fits.fz) the ones of the list of this night and telescope
raw_files() {
  local night="${DATE:0:4}/${DATE:4:2}/${DATE:6:2}"
  if [[ -n "${RAW_LIST:-}" ]]; then
    # lines without comments and blanks, relative paths completed, the files of the night and telescope
    local file
    sed -e 's/#.*//' -e 's/[[:space:]]//g' "${RAW_LIST}" \
      | awk -v base="${RAW_DIRECTORY}/cs$1/" 'NF {print ($0 ~ /^\// ? $0 : base $0)}' \
      | { grep -E "/${night}/(.*/)?SST1M$1_[^/]*\.fits\.fz$" || true; } | sort -u \
      | while IFS= read -r file; do
          if [[ -f "${file}" ]]; then
            echo "${file}"
          else
            echo "warning: ${file} (RAW_LIST) does not exist" >&2
          fi
        done
    return 0
  fi
  local directory="${RAW_DIRECTORY}/cs$1/data/raw/${night}"
  [[ -d "${directory}" ]] || return 0
  find "${directory}" -name "SST1M$1_*.fits.fz" | sort
}

# nights (YYYYMMDD) of the raw files of the list $1 (directories .../YYYY/MM/DD/), sorted
list_nights() {
  grep -oE '/[0-9]{4}/[0-9]{2}/[0-9]{2}/' "$1" | tr -d / | sort -u
}
