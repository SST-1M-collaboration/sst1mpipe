# Daily analysis of the SST-1M data

The data of a night, R0 -> DL1 -> DL1 stereo -> DL2 -> DL3, with the models and IRFs of the simulations
(`scripts/camk/submit_mc.bash`). It replaces `daily_analysis.py`.

## On Calculus (one machine, 32 cores): `run_daily.bash`

```bash
bash scripts/daily/run_daily.bash                # the night before (yesterday)
DATE=20260120 bash scripts/daily/run_daily.bash  # a given night
```

- **R0 -> DL1:** `N_PARALLEL` raw files at a time (default: the cores of the machine, `nproc`).
- **DL1 stereo:** run when both telescopes observed.
- **DL2 and DL3:** `N_PARALLEL` files at a time.
- **Outputs:** in `DAILY_DIRECTORY/<DATE>/` (default `/data/work/analysis/Daily_analysis/daily_data/<DATE>/`).
- **Logs:** one per step and per raw file, in `<DATE>/logs/`.
- **Failures:** a file that fails is reported and the next steps go on with the files that succeeded.
- **Empty night:** a night without raw files exits without error.
- **One run at a time:** a lock file prevents two runs from overlapping.

Every morning with cron (`crontab`, at 10:17): set `SST1MPIPE_ROOT` in the file, then

```bash
(crontab -l 2>/dev/null; grep -v '^#' scripts/daily/crontab) | crontab -
```

The log of the cron runs is `~/logs/sst1mpipe-daily.log`. The conda environment (`CONDA_ENV`, default
`sst1m-dev`) is activated by the scripts. Conda is looked for in the `PATH`, then in `CONDA_ROOT`, then
in `~/miniforge3`, `~/miniconda3` and `~/anaconda3`.

Settings are in `daily.env`; every variable can also be set in the environment:

| Variable | Default (Calculus) |
|---|---|
| `RAW_DIRECTORY` | `/net`: raw files in `cs<tel>/data/raw/YYYY/MM/DD/` |
| `DAILY_DIRECTORY` | `/data/work/analysis/Daily_analysis/daily_data` |
| `MC_DIRECTORY` | `/data/work/analysis/MC/crab_paper_NSBmed4` |
| `MODELS_NAME` | `ze20_30` |
| `IRF_CUTS` | `optimized` |
| `MODELS_DIRECTORY` | `MC_DIRECTORY/models/MODELS_NAME`, with `mono/` and `stereo/` |
| `DATA_IRF_DIRECTORY` | `MC_DIRECTORY/irf/MODELS_NAME/IRF_CUTS`, with `tel1/`, `tel2/` and `stereo/` |
| `N_PARALLEL` | `nproc` |
| `DL2_N_JOBS` | 1 thread of the random forests per DL2 file |

Copy the models and IRFs from camk into `MC_DIRECTORY` (same layout), e.g. on Calculus:

```bash
rsync -a <user>@ui.cta.camk.edu.pl:/work/se2/<user>/sst1mpipe/mc/crab_paper_NSBmed4/{models,irf} /data/work/analysis/MC/crab_paper_NSBmed4/
```

## On camk (SLURM): `submit_data.bash`

```bash
DATE=20260120 bash scripts/daily/submit_data.bash
```

The same steps run as SLURM jobs: one array task per raw file for R0 -> DL1, then one job per step. The
settings are in `camk.env` (the camk paths, then `daily.env`), and the logs in `~/logs`.

## Steps

| Script | |
|---|---|
| `data_r0_dl1.bash` | `sst1mpipe-monitoring`, `sst1mpipe-process` of one raw file (argument, or the task of the SLURM job array) |
| `data_dl1_stereo.bash` | `ctapipe-merge` of the night of each telescope, `sst1mpipe-stereo` (HillasReconstructor) |
| `data_dl1_dl2.bash` | `sst1mpipe-apply-models` on every DL1 file: cs1 and cs2 with the mono models, stereo with the stereo models |
| `data_dl2_dl3.bash` | `sst1mpipe-dl2-dl3` on every DL2 file with the IRFs of `DATA_IRF_DIRECTORY` (cs1: `tel1`, cs2: `tel2`, stereo: `stereo`, closest zenith angle), then the gammapy index files of `dl3/<cs1, cs2, stereo>/` |

DL3 with gammapy: `DataStore.from_dir("<DATE>/dl3/stereo").get_observations(required_irf="point-like")`.
TODO: dead time and GTI, issue #252.
