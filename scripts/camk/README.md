# SST-1M pipeline on camk (SLURM)

Paths and settings: `sst1mpipe.env` (default simulations: `/work/se3/jurysek/sst1m_crab_paper_mc`, `<particle>/zenith_<zenith>deg/`) (every variable can be set in the environment instead). The
submit scripts run on the login node and submit the jobs with their dependencies; logs in `~/logs`.

```bash
bash scripts/camk/submit_mc.bash                  # simulations R0 -> IRFs
DATE=20260120 bash scripts/camk/submit_data.bash  # data of a night R0 -> DL3
```

| Simulations | |
|---|---|
| `mc_r0_dl1.bash` | `sst1mpipe-process`, `MC_FILES_PER_TASK` simtel files per array task, per particle and zenith angle (`MC_ZENITHS`) |
| `mc_merge.bash` | `ctapipe-merge`: training (`MC_TRAIN_FRACTION`) and test files of each zenith angle in `merged/ze<z>/`; with several training zenith angles (`MC_TRAIN_ZENITHS`), their training files in `merged/train_<MODELS_NAME>/` |
| `mc_train.bash` | `sst1mpipe-train-*`: energy, particle classifier, disp in `models/<MODELS_NAME>/` (the NSB level, feature `cos_zenith`), with plots |
| `mc_dl1_dl2.bash` | `sst1mpipe-apply-models` on the test files of each zenith angle: `dl2/<MODELS_NAME>/ze<z>/stereo` (both telescopes), `tel1`, `tel2` (one telescope alone); the random forests use the `MC_DL2_CPUS` CPUs of each task (default 16, `--n-jobs`) |
| `mc_irf.bash` | `ctapipe-optimize-event-selection`, `ctapipe-compute-irf`: point-like and full enclosure IRFs per zenith angle and mode (`tel1`, `tel2`, `stereo`, `strict_stereo`), in `irf/<MODELS_NAME>/<IRF_CUTS>/<mode>/ze<z>/`, with their figures (`plots/`, `sst1mpipe-plot-irfs`) |
| `mc_irf_comparison.bash` | figures of the four modes together, in `irf/<MODELS_NAME>/<IRF_CUTS>/comparison/ze<z>/` |

`MODELS_NAME`: the zenith angles of the training, `ze20` or `ze20_30` (`MC_TRAIN_ZENITHS`, default
`MC_ZENITHS`): the models, their DL2 and IRFs of different trainings do not overwrite each other. E.g.
a new zenith angle (ze30) and models trained on ze20 and ze30, the files of ze20 being merged:

```bash
MC_ZENITHS="30" MC_TRAIN_ZENITHS="20 30" bash scripts/camk/submit_mc.bash
```

and the DL2 and IRFs of ze20 with these models: `MC_START=dl2 MC_ZENITHS="20" MC_TRAIN_ZENITHS="20 30"`.
The simulations of the zenith angles must have different CORSIKA runs (obs_id): `mc_merge.bash` stops otherwise.

| Data | |
|---|---|
| `data_r0_dl1.bash` | `sst1mpipe-monitoring`, `sst1mpipe-process`, one raw file per array task, per telescope (raw files: `RAW_DIRECTORY`, default `/work/se1/sst1m_data/cs<tel>/data/raw/YYYY/MM/DD`) |
| `data_dl1_stereo.bash` | `ctapipe-merge` of the night of each telescope, `sst1mpipe-stereo` |
| `data_dl1_dl2.bash` | `sst1mpipe-apply-models` on every DL1 file (cs1, cs2, stereo) |
| `data_dl2_dl3.bash` | `sst1mpipe-dl2-dl3` on every DL2 file with the IRFs of `DATA_IRF_DIRECTORY` (cs1: `tel1`, cs2: `tel2`, stereo: `strict_stereo`, closest zenith angle), index files of gammapy per analysis in `dl3/<cs1, cs2, stereo>/` |


IRFs: `IRF_CUTS=optimized` (best sensitivity) or `paper` (fixed cuts of arXiv:2506.01733: 60 % gamma
efficiency, 80 % theta containment); background weighted to protons + helium (`IRF_BACKGROUND=p_he`,
or `protons`). `MC_START=dl2` / `irf`: from the DL2 / IRFs of a previous submission, e.g.

```bash
MC_START=irf IRF_CUTS=paper bash scripts/camk/submit_mc.bash
```

DL3 (gammapy): `DataStore.from_dir("$DATA_DIRECTORY/dl3/stereo").get_observations(required_irf="point-like")`. TODO: dead time and GTI, issue #252.

Comparison with the published performance (`ref:`: IRF directories of `sst1mpipe/resources/reference`,
e.g. the SST-1M Crab paper, `crab_paper_2025/<mode>/ze<z>`, see `sst1mpipe-plot-irfs --help`):

```bash
sst1mpipe-plot-irfs --irf strict_stereo=$MC_DIRECTORY/irf/ze20/paper/strict_stereo/ze20 \
    --reference "paper 30°=ref:crab_paper_2025/strict_stereo/ze30" \
    --reference "paper 20°=ref:crab_paper_2025/strict_stereo/ze20" -o plots_paper
```
