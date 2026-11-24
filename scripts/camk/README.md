# SST-1M pipeline on camk (SLURM)

Paths and settings: `sst1mpipe.env` (default simulations: `/work/se3/jurysek/sst1m_crab_paper_mc`, `<particle>/zenith_<zenith>deg/`) (every variable can be set in the environment instead). The
submit scripts run on the login node and submit the jobs with their dependencies; logs in `~/logs`.

```bash
bash scripts/camk/submit_mc.bash  # simulations R0 -> IRFs
```

The data of a night (R0 -> DL3) are processed by the scripts of `scripts/daily` (camk or Calculus).

| Simulations | |
|---|---|
| `mc_r0_dl1.bash` | `sst1mpipe-process`, `MC_FILES_PER_TASK` simtel files per array task, per particle and zenith angle (`MC_ZENITHS`) |
| `mc_merge.bash` | `ctapipe-merge`, one array task per merged file (5 per zenith angle): training (`MC_TRAIN_FRACTION`) and test files of each zenith angle in `merged/ze<z>/` |
| `mc_merge_train.bash` | with several training zenith angles (`MC_TRAIN_ZENITHS`): their training files in `merged/train_<MODELS_NAME>/`, one array task per particle |
| `mc_train.bash` | `sst1mpipe-train-*`: energy, particle classifier, disp in `models/<MODELS_NAME>/mono/` (image parameters) and `models/<MODELS_NAME>/stereo/` (and the impact distance and `h_max` of the HillasReconstructor, events with both telescopes), with plots |
| `mc_dl1_dl2.bash` | `sst1mpipe-apply-models` on the test files of each zenith angle: `dl2/<MODELS_NAME>/ze<z>/stereo` (stereo models, events with both telescopes), `tel1`, `tel2` (mono models, one telescope alone); the random forests use the `MC_DL2_CPUS` CPUs of each task (default 16, `--n-jobs`) |
| `mc_irf.bash` | `sst1mpipe-optimize-event-selection`, `sst1mpipe-compute-irf` (ctapipe IRF tools with the Crab of MAGIC and protons + helium of the paper): point-like and full enclosure IRFs per zenith angle and mode (`tel1`, `tel2`, `stereo`), in `irf/<MODELS_NAME>/<IRF_CUTS>/<mode>/ze<z>/`, with their figures (`plots/`, `sst1mpipe-plot-irfs`) |
| `mc_irf_comparison.bash` | figures of the three modes together, in `irf/<MODELS_NAME>/<IRF_CUTS>/comparison/ze<z>/` |

`MODELS_NAME`: the zenith angles of the training, `ze20` or `ze20_30` (`MC_TRAIN_ZENITHS`, default
`MC_ZENITHS`): the models, their DL2 and IRFs of different trainings do not overwrite each other. E.g.
a new zenith angle (ze30) and models trained on ze20 and ze30, the files of ze20 being merged:

```bash
MC_ZENITHS="30" MC_TRAIN_ZENITHS="20 30" bash scripts/camk/submit_mc.bash
```

and the DL2 and IRFs of ze20 with these models: `MC_START=dl2 MC_ZENITHS="20" MC_TRAIN_ZENITHS="20 30"`.
The simulations of the zenith angles have the same CORSIKA runs (obs_id): in the training files of
several zenith angles, the obs_ids of each zenith angle get the offset zenith × 10000000 (`offset_obs_ids.py`).

IRFs: `IRF_CUTS=optimized` (best sensitivity) or `paper` (fixed cuts of arXiv:2506.01733: 60 % gamma
efficiency per energy bin, one theta cut of 80 % containment); background weighted to protons + helium (`IRF_BACKGROUND=p_he`,
or `protons`). `MC_START=dl2` / `irf`: from the DL2 / IRFs of a previous submission, e.g.

```bash
MC_START=irf IRF_CUTS=paper bash scripts/camk/submit_mc.bash
```

Comparison with the published performance (`ref:`: IRF directories of `sst1mpipe/resources/reference`,
e.g. the SST-1M Crab paper, `crab_paper_2025/<mode>/ze<z>`, see `sst1mpipe-plot-irfs --help`):

```bash
sst1mpipe-plot-irfs --irf stereo=$MC_DIRECTORY/irf/ze20/paper/stereo/ze20 \
    --reference "paper 30°=ref:crab_paper_2025/stereo/ze30" \
    --reference "paper 20°=ref:crab_paper_2025/stereo/ze20" -o plots_paper
```
