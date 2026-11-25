#!/bin/bash

##########################################################################
#  Stereo direction: the combiners compared on the same models (camk)    #
##########################################################################
# The stereo models of a previous submission (models/<MODELS_NAME>/stereo) applied to the test files with
# four stereo combinations of the direction (sst1mpipe-apply-models --ApplyModels.disp_stereo_combiner):
#   stereo_mars        StereoDispSignCrossingCombiner (former StereoMARSCombiner): signs of the crossing
#                      of the major axes, events with at least one valid telescope (as the IRFs before)
#   stereo_mars_2tel   the same, events with both telescopes valid
#   stereo_paper_1tel  SST-1M Crab paper: closest combination of the signs, events with at least one telescope
#   stereo_paper       SST-1M Crab paper: closest combination of the signs, both telescopes (StereoDispSignClosestCombiner)
# then their IRFs (IRF_CUTS, default paper) and the comparison:
#   irf/<MODELS_NAME>/<IRF_CUTS>/combiners/ze<zenith>/: figures of the four versions with the paper, and
#   stereo_direction_<zenith>.txt/png: theta68 without cuts, per number of telescopes and sign method
#
#   MC_ZENITHS="30" MC_TRAIN_ZENITHS="30" bash scripts/camk/submit_stereo_combiners.bash
# MC_START: train (default: the disp models, mono and stereo, are retrained first: the models trained
# with the former StereoMARSCombiner cannot be read any more), dl2 (models already retrained) or irf (DL2 of a previous
# submission). STEREO_VERSIONS: other versions (subset of the four).

export SST1MPIPE_ROOT="${SST1MPIPE_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
export IRF_CUTS="${IRF_CUTS:-paper}"
CONDA_ENV="" source "${SST1MPIPE_ROOT}/scripts/camk/common.bash"
SCRIPTS="${SST1MPIPE_ROOT}/scripts/camk"

export STEREO_VERSIONS="${STEREO_VERSIONS:-stereo_mars stereo_mars_2tel stereo_paper_1tel stereo_paper}"
export MC_DL2_VERSIONS="${STEREO_VERSIONS}" IRF_MODES="${STEREO_VERSIONS}"
read -ra VERSIONS <<< "${STEREO_VERSIONS}"
read -ra ZENITHS <<< "${MC_ZENITHS}"
N=$(( ${#VERSIONS[@]} * ${#ZENITHS[@]} ))
if [[ "${MC_START:-train}" == train ]]; then
  [[ -f "$(train_directory)/gamma_diffuse_train.dl1.h5" ]] || { echo "No training files in $(train_directory)" >&2; exit 1; }
elif ! compgen -G "${MC_DIRECTORY}/models/${MODELS_NAME}/stereo/*.pkl" > /dev/null; then
  echo "No stereo models in ${MC_DIRECTORY}/models/${MODELS_NAME}/stereo (MODELS_NAME from MC_TRAIN_ZENITHS=\"${MC_TRAIN_ZENITHS}\")" >&2
  exit 1
fi

MC_START="${MC_START:-train}"
[[ " train dl2 irf " == *" ${MC_START} "* ]] || { echo "MC_START: train, dl2 or irf" >&2; exit 1; }
DEPENDENCY=()
if [[ "${MC_START}" == train ]]; then
  # disp reconstructors (tasks 2: mono, 5: stereo of mc_train.bash)
  TRAIN=$(sbatch --parsable --array=2,5 "${SCRIPTS}/mc_train.bash")
  DEPENDENCY=(--dependency="afterok:${TRAIN}")
  echo "training of the disp models: job ${TRAIN}"
fi
if [[ "${MC_START}" != irf ]]; then
  # 3 test files per version and zenith angle
  DL2=$(sbatch --parsable "${DEPENDENCY[@]}" --array="0-$((3 * N - 1))" --cpus-per-task="${MC_DL2_CPUS}" \
    --mem-per-cpu="${MC_DL2_MEM_PER_CPU}" "${SCRIPTS}/mc_dl1_dl2.bash")
  DEPENDENCY=(--dependency="afterok:${DL2}")
  echo "DL2 (${STEREO_VERSIONS}): job ${DL2}"
fi
IRF=$(sbatch --parsable "${DEPENDENCY[@]}" --array="0-$((N - 1))" "${SCRIPTS}/mc_irf.bash")
COMPARISON=$(sbatch --parsable --dependency="afterany:${IRF}" --array="0-$(( ${#ZENITHS[@]} - 1 ))" \
  "${SCRIPTS}/mc_stereo_combiners_comparison.bash")
echo "IRFs (${IRF_CUTS} cuts): job ${IRF}, comparison ${COMPARISON}:"
echo "  ${MC_DIRECTORY}/irf/${MODELS_NAME}/${IRF_CUTS}/combiners/ze<zenith>/"
