#!/bin/bash

################################
#  DL1 -> DL2 of simulations   #
################################
# Test files (point gammas, diffuse gammas, protons) of each zenith angle of MC_ZENITHS, with the
# models of models/<MODELS_NAME>, in three versions: dl2/<MODELS_NAME>/ze<zenith>/<version>/
# - stereo: the stereo models (models/<MODELS_NAME>/stereo): the events with both telescopes, combined
# - tel1, tel2: the mono models (models/<MODELS_NAME>/mono), one telescope alone (its mono reconstruction)
# Task 3 n i + 3 j + k: zenith angle i, version j (n versions, MC_DL2_VERSIONS), file k.
# The random forests use the CPUs of the task (--n-jobs): MC_DL2_CPUS (submit_mc.bash, default 16).

#SBATCH --job-name=mc-dl1-dl2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-08:00:00
#SBATCH --qos=long
#SBATCH --partition=pctacamk
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/home/%u/logs/%x-%A_%a.out
#SBATCH --error=/home/%u/logs/%x-%A_%a.err

source "${SST1MPIPE_ROOT:-${SLURM_SUBMIT_DIR:-.}}/scripts/camk/common.bash"

read -ra ZENITHS <<< "${MC_ZENITHS}"
NAMES=(gamma_point gamma_diffuse_test proton_test)
# MC_DL2_VERSIONS: the versions (default stereo tel1 tel2); stereo_<combiner>: the stereo models with
# another stereo combiner of the direction (submit_stereo_combiners.bash)
read -ra VERSIONS <<< "${MC_DL2_VERSIONS:-stereo tel1 tel2}"
N_VERSIONS=${#VERSIONS[@]}
TASK=${SLURM_ARRAY_TASK_ID:-0}
ZENITH=${ZENITHS[TASK / (3 * N_VERSIONS)]}
VERSION=${VERSIONS[TASK % (3 * N_VERSIONS) / 3]}
NAME=${NAMES[TASK % 3]}
# telescopes of the models, stereo combiner of the direction
COMBINER=()
case "${VERSION}" in
  tel1) TELESCOPES="[1]" ;;
  tel2) TELESCOPES="[2]" ;;
  stereo) TELESCOPES="[]" ;;
  stereo_*)
    TELESCOPES="[]"
    case "${VERSION#stereo_}" in
      # StereoDispSignCrossingCombiner (former StereoMARSCombiner): signs of the crossing of the major axes
      mars) COMBINER=(StereoDispSignCrossingCombiner 1) ;;
      mars_2tel) COMBINER=(StereoDispSignCrossingCombiner 2) ;;
      # StereoDispSignClosestCombiner (SST-1M Crab paper): closest combination of the signs
      paper_1tel) COMBINER=(StereoDispSignClosestCombiner 1) ;;
      paper) COMBINER=(StereoDispSignClosestCombiner 2) ;;
      *) echo "Unknown version ${VERSION}: stereo_mars, stereo_mars_2tel, stereo_paper_1tel or stereo_paper" >&2; exit 1 ;;
    esac
    # minimum number of valid telescopes of the stereo events
    COMBINER=(--ApplyModels.disp_stereo_combiner="${COMBINER[0]}" --"${COMBINER[0]}".weights=intensity
              --"${COMBINER[0]}".min_telescopes="${COMBINER[1]}")
    ;;
  *) echo "Unknown version ${VERSION}" >&2; exit 1 ;;
esac
OUTPUT="${MC_DIRECTORY}/dl2/${MODELS_NAME}/ze${ZENITH}/${VERSION}/${NAME}.dl2.h5"
mkdir -p "$(dirname "${OUTPUT}")"
SET=mono
[[ "${VERSION}" == stereo* ]] && SET=stereo
mapfile -t MODELS < <(models_args "${MC_DIRECTORY}/models/${MODELS_NAME}/${SET}")
(( ${#MODELS[@]} > 0 )) || exit 1
sst1mpipe-apply-models -i "${MC_DIRECTORY}/merged/ze${ZENITH}/${NAME}.dl1.h5" "${MODELS[@]}" -o "${OUTPUT}" \
  --ApplyModels.telescopes="${TELESCOPES}" "${COMBINER[@]}" --n-jobs "${SLURM_CPUS_PER_TASK:-1}" --overwrite --no-progress
