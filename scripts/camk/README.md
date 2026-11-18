# SST-1M pipeline on camk (SLURM)

Paths and settings: `sst1mpipe.env` (every variable can be set in the environment instead). The
submit scripts run on the login node and submit the jobs with their dependencies; logs in `~/logs`.

```bash
bash scripts/camk/submit_mc.bash                  # simulations R0 -> IRFs
DATE=20260120 bash scripts/camk/submit_data.bash  # data of a night R0 -> DL2
```

| Simulations | |
|---|---|
| `mc_r0_dl1.bash` | `sst1mpipe-process`, one simtel file per array task, per particle |
| `mc_merge.bash` | `ctapipe-merge`: point gammas, diffuse gammas and protons split in training / test files (`MC_TRAIN_FRACTION`) |
| `mc_train.bash` | `sst1mpipe-train-*`: energy, particle classifier, disp in `models/<MC_MODEL_NODE>`, with plots |
| `mc_dl1_dl2.bash` | `sst1mpipe-apply-models` on the test files |
| `mc_irf.bash` | `ctapipe-optimize-event-selection`, `ctapipe-compute-irf`: point-like and full enclosure IRFs, in `IRF_CONDA_ENV` |

| Data | |
|---|---|
| `data_r0_dl1.bash` | `sst1mpipe-monitoring`, `sst1mpipe-process`, one raw file per array task, per telescope |
| `data_dl1_stereo.bash` | `sst1mpipe_data_dl1_dl1_stereo` |
| `data_dl1_dl2.bash` | `sst1mpipe-apply-models` on every DL1 file (cs1, cs2, stereo) |

`IRF_CONDA_ENV`: the IRF tools of ctapipe need pyirf ~0.14 (astropy >= 7), sst1mpipe pyirf ~0.8
(astropy < 7). TODO: DL3 of the data.
