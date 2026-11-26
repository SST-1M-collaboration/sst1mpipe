# Processing of a dataset

A dataset is the list of the runs of an analysis in an ECSV file (`sst1mpipe.io.dataset`), e.g.
`sst1mpipe/resources/datasets/crab_paper_2025.ecsv` (git LFS: `git lfs pull -I 'sst1mpipe/resources/datasets/*'`).
`submit_dataset.bash` processes its runs on camk (SLURM) by reusing the other scripts:

- **each night**: `scripts/daily/submit_data.bash` with the raw files of the dataset of this night
  (`RAW_LIST`): monitoring and R0 -> DL1 per raw file, DL1 stereo, DL2, DL3;
- **models and IRFs**: those of the simulations, `scripts/camk/submit_mc.bash` (`MC_DIRECTORY`,
  `MODELS_NAME`, `IRF_CUTS` of `scripts/daily/camk.env` / `daily.env`);
- **all the nights**: `dataset_dl3_index.bash`, the index files of gammapy (`create_hdu_indexes`) of the
  DL3 files of all the nights.

```bash
DATASET=crab_paper_2025 DATA_STOP=dl1 bash scripts/dataset/submit_dataset.bash   # monitoring, R0 -> DL1
DATASET=crab_paper_2025 bash scripts/dataset/submit_dataset.bash                 # up to DL3
DATASET=crab_paper_2025 DATASET_NIGHTS="20240110 20240111" bash scripts/dataset/submit_dataset.bash
DATASET=my_runs.ecsv MODELS_NAME=ze20_30 IRF_CUTS=paper bash scripts/dataset/submit_dataset.bash
```

| Variable | |
|---|---|
| `DATASET` | name of a dataset of sst1mpipe, or ECSV file |
| `DATASET_ROOT` | root of the raw files: the `root` of the dataset if it exists, else `RAW_DIRECTORY` (camk: `/work/se1/sst1m_data`) |
| `DATASET_NIGHTS` | only these nights (`YYYYMMDD`, separated by spaces) |
| `DATA_STOP` | last step: `dl1`, `stereo`, `dl2` or `dl3` (default) |
| `DATASETS_DIRECTORY` | outputs: `<DATASETS_DIRECTORY>/<name>/` (default: `/work/se2/<user>/sst1mpipe/datasets`) |

Outputs in `<DATASETS_DIRECTORY>/<name>/`:

- `dataset/`: the dataset, the list of the raw files, the configurations, `settings.txt` (settings,
  sst1mpipe version and branch, nights), to reproduce the analysis;
- `<night>/`: the outputs of the night, as the daily analysis (`lists`, `monitoring`, `dl1`, `dl2`, `dl3`);
- `dl3/<cs1, cs2, stereo>/`: index files of all the nights, for gammapy:
  `DataStore.from_dir("<DATASETS_DIRECTORY>/<name>/dl3/stereo")`.

The steps of a night that fail are skipped by the next ones (`afterany`); logs in `~/logs`.
