"""
Stereo direction: where is the angular resolution lost compared with the SST-1M Crab paper?

Reads a stereo DL2 file of point gammas (sst1mpipe-apply-models with the stereo models) and computes
the 68 % containment of the angular distance to the true direction, per true energy bin (5 per decade,
as the paper), for:

  pipeline       the stereo direction of the pipeline (disp, StereoMARSCombiner) of all the valid events,
                 as in the IRFs: events with one valid telescope included
  1 telescope    the events of the pipeline with one valid telescope only (the other one failed the quality
                 cuts of the models: intensity > 45 p.e., leakage2 < 0.7)
  2 telescopes   the events of the pipeline with both telescopes valid (paper and sst1mpipe < 0.9: stereo
                 events are the events with both telescopes passing the cuts)
  paper signs    events with 2 telescopes, signs of the paper / sst1mpipe < 0.9 (MARS): of the 4 sign
                 combinations, the one whose two positions are the closest; intensity-weighted mean
  true signs     events with 2 telescopes, the sign combination closest to the true direction: the best
                 that the disp norm allows (limit of any sign choice)
  Hillas         events with 2 telescopes, the HillasReconstructor direction (crossing of the major axes)

No gammaness or theta cut (pure reconstruction). Also writes the fraction of events with one telescope
and the fraction of events whose sign combination differs from the true one, per method.

    python stereo_direction_diagnostic.py <dl2 dir>/stereo/gamma_point.dl2.h5 <output dir> [label]
"""
import sys
from pathlib import Path

import astropy.units as u
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import join, vstack
from ctapipe.io import read_table
from ctapipe.reco.preprocessing import telescope_to_horizontal

path, out = Path(sys.argv[1]), Path(sys.argv[2])
label = sys.argv[3] if len(sys.argv) > 3 else path.parent.name
out.mkdir(parents=True, exist_ok=True)
KEY = ["obs_id", "event_id"]


def tel_ids(path):
    import tables
    with tables.open_file(path) as f:
        return sorted(int(n._v_name.split("_")[1]) for n in f.get_node("/dl2/event/telescope/disp/disp"))


def unit_vector(alt, az):
    alt, az = np.deg2rad(alt), np.deg2rad(az)
    return np.stack([np.cos(alt) * np.cos(az), np.cos(alt) * np.sin(az), np.sin(alt)], axis=-1)


def separation(alt1, az1, alt2, az2):
    dot = np.sum(unit_vector(alt1, az1) * unit_vector(alt2, az2), axis=-1)
    return np.rad2deg(np.arccos(np.clip(dot, -1, 1)))


def mean_direction(alts, azs, weights):
    v = np.sum(unit_vector(alts, azs) * weights[..., None], axis=0)
    v /= np.linalg.norm(v, axis=-1, keepdims=True)
    return np.rad2deg(np.arcsin(v[..., 2])), np.rad2deg(np.arctan2(v[..., 1], v[..., 0])) % 360


# ---------------------------------------------------------------- per telescope
tels = tel_ids(path)
per_tel = {}
for tel in tels:
    name = f"tel_{tel:03d}"
    t = read_table(path, f"/dl1/event/telescope/parameters/{name}")[KEY + ["hillas_fov_lon", "hillas_fov_lat", "hillas_psi", "hillas_intensity"]]
    t = join(t, read_table(path, f"/dl2/event/telescope/disp/disp/{name}")[KEY + ["disp_tel_parameter"]], keys=KEY)
    t = join(t, read_table(path, f"/dl2/event/telescope/geometry/disp/{name}")[KEY + ["disp_tel_is_valid"]], keys=KEY)
    pointing = read_table(path, f"/configuration/telescope/pointing/{name}")
    alt0 = u.Quantity(pointing["telescope_pointing_altitude"][0], u.rad).to(u.deg)
    az0 = u.Quantity(pointing["telescope_pointing_azimuth"][0], u.rad).to(u.deg)
    lon, lat = t["hillas_fov_lon"].quantity.to(u.deg), t["hillas_fov_lat"].quantity.to(u.deg)
    psi, norm = t["hillas_psi"].quantity.to(u.rad), np.abs(t["disp_tel_parameter"].quantity.to(u.deg))
    for sign, s in (("p", 1), ("m", -1)):
        alt, az = telescope_to_horizontal(lon=lon + s * norm * np.cos(psi), lat=lat + s * norm * np.sin(psi),
                                          pointing_alt=alt0, pointing_az=az0)
        t[f"alt_{sign}"], t[f"az_{sign}"] = alt.to_value(u.deg), az.to_value(u.deg)
    t["predicted"] = np.where(np.asarray(t["disp_tel_parameter"]) < 0, "m", "p")
    keep = ["obs_id", "event_id", "hillas_intensity", "disp_tel_is_valid", "alt_p", "az_p", "alt_m", "az_m", "predicted"]
    per_tel[tel] = t[keep]
    for c in keep[2:]:
        per_tel[tel].rename_column(c, f"{c}_{tel}")

# ---------------------------------------------------------------- per event
events = read_table(path, "/dl2/event/subarray/geometry/disp")[KEY + ["disp_alt", "disp_az", "disp_is_valid"]]
events = join(events, read_table(path, "/dl2/event/subarray/geometry/HillasReconstructor")[KEY + ["HillasReconstructor_alt", "HillasReconstructor_az", "HillasReconstructor_is_valid"]], keys=KEY)
events = join(events, read_table(path, "/simulation/event/subarray/shower")[KEY + ["true_energy", "true_alt", "true_az"]], keys=KEY)
for tel in tels:
    events = join(events, per_tel[tel], keys=KEY, join_type="left")
events = events.filled(np.nan) if events.masked else events
for tel in tels:
    v = events[f"disp_tel_is_valid_{tel}"]
    events[f"disp_tel_is_valid_{tel}"] = np.nan_to_num(np.asarray(v, dtype=float), nan=0).astype(bool)

energy = events["true_energy"].quantity.to_value(u.TeV)
true_alt, true_az = events["true_alt"].quantity.to_value(u.deg), events["true_az"].quantity.to_value(u.deg)
n_valid = sum(np.asarray(events[f"disp_tel_is_valid_{tel}"]) for tel in tels)
valid = np.asarray(events["disp_is_valid"], dtype=bool)
two = valid & (n_valid >= 2)
one = valid & (n_valid == 1)

theta = {"pipeline": separation(events["disp_alt"], events["disp_az"], true_alt, true_az)}
masks = {"pipeline": valid, "1 telescope": one, "2 telescopes": two}
theta["1 telescope"] = theta["2 telescopes"] = theta["pipeline"]

if len(tels) == 2:
    a, b = tels
    combos = ["pp", "mm", "pm", "mp"]
    pos = {c: ((events[f"alt_{c[0]}_{a}"], events[f"az_{c[0]}_{a}"]), (events[f"alt_{c[1]}_{b}"], events[f"az_{c[1]}_{b}"])) for c in combos}
    w = np.stack([np.asarray(events[f"hillas_intensity_{a}"], float), np.asarray(events[f"hillas_intensity_{b}"], float)])

    def direction(choice):
        alts = np.stack([np.choose(choice, [np.asarray(pos[c][i][0], float) for c in combos]) for i in (0, 1)])
        azs = np.stack([np.choose(choice, [np.asarray(pos[c][i][1], float) for c in combos]) for i in (0, 1)])
        return mean_direction(alts, azs, w)

    distance = np.stack([separation(*pos[c][0], *pos[c][1]) for c in combos])
    paper_choice = np.argmin(np.nan_to_num(distance, nan=1e9), axis=0)
    true_choice = np.argmin(np.nan_to_num(np.stack([separation(*direction(np.full(len(events), i)), true_alt, true_az)
                                                    for i in range(4)]), nan=1e9), axis=0)
    theta["paper signs"] = separation(*direction(paper_choice), true_alt, true_az)
    theta["true signs"] = separation(*direction(true_choice), true_alt, true_az)
    masks["paper signs"] = masks["true signs"] = two
    # sign combination of the pipeline: the closest combination to its direction
    pipeline_choice = np.argmin(np.nan_to_num(np.stack([separation(*direction(np.full(len(events), i)), events["disp_alt"], events["disp_az"])
                                                        for i in range(4)]), nan=1e9), axis=0)
    predicted_choice = np.array([combos.index(x + y) if isinstance(x, str) and isinstance(y, str) and (x + y) in combos else 0
                                 for x, y in zip(events[f"predicted_{a}"], events[f"predicted_{b}"])])
theta["Hillas"] = separation(events["HillasReconstructor_alt"], events["HillasReconstructor_az"], true_alt, true_az)
masks["Hillas"] = two & np.asarray(events["HillasReconstructor_is_valid"], dtype=bool)

# ---------------------------------------------------------------- per energy
bins = np.logspace(-1, 3, 21)
centers = np.sqrt(bins[:-1] * bins[1:])
index = np.digitize(energy, bins) - 1


def per_bin(values, mask, q=68, min_events=10):
    r = np.full(len(centers), np.nan)
    for i in range(len(centers)):
        v = values[mask & (index == i)]
        v = v[np.isfinite(v)]
        if len(v) >= min_events:
            r[i] = np.percentile(v, q)
    return r


lines = [f"# {path}", "# theta68 (deg) per true energy bin, no gammaness/theta cut; n: events", ""]
header = f"{'E (TeV)':>8} " + " ".join(f"{k:>13}" for k in theta) + f" {'1-tel frac':>10}"
if len(tels) == 2:
    header += f" {'wrong sign: pipeline':>21} {'paper':>6} {'predicted':>9}"
lines.append(header)
results = {k: per_bin(theta[k], masks[k]) for k in theta}
for i, e in enumerate(centers):
    in_bin = index == i
    n = np.count_nonzero(valid & in_bin)
    if n < 10:
        continue
    row = f"{e:8.3g} " + " ".join(f"{results[k][i]:13.4f}" for k in theta) + f" {np.count_nonzero(one & in_bin) / n:10.2f}"
    if len(tels) == 2:
        m = two & in_bin
        nm = max(np.count_nonzero(m), 1)
        row += (f" {np.count_nonzero(pipeline_choice[m] != true_choice[m]) / nm:21.2f}"
                f" {np.count_nonzero(paper_choice[m] != true_choice[m]) / nm:6.2f}"
                f" {np.count_nonzero(predicted_choice[m] != true_choice[m]) / nm:9.2f}")
    lines.append(row)
lines.append("")
lines.append("events: " + ", ".join(f"{k} {np.count_nonzero(masks[k])}" for k in masks))
text = "\n".join(lines)
(out / f"stereo_direction_{label}.txt").write_text(text)
print(text)

# reference: the paper, stereo 30 deg
try:
    from sst1mpipe.performance.irf_plots import REFERENCE_DIRECTORY, read_irf_directory
    ref = read_irf_directory(REFERENCE_DIRECTORY / "crab_paper_2025" / "stereo" / "ze30")["angular"]
    ref_lo, ref_hi = (np.ravel(u.Quantity(ref[c], u.TeV).value) for c in ("ENERG_LO", "ENERG_HI"))
    ref_v = np.ravel(ref["ANGULAR_RESOLUTION_68"].to_value(u.deg))
except Exception:  # noqa: BLE001
    ref = None

fig, ax = plt.subplots(figsize=(7, 4.5))
styles = {"pipeline": ("k", "o"), "1 telescope": ("#eb6834", "v"), "2 telescopes": ("#2a78d6", "o"),
          "paper signs": ("#1baf7a", "s"), "true signs": ("#6250d6", "^"), "Hillas": ("#eda100", "D")}
for k, v in results.items():
    color, marker = styles[k]
    ax.plot(centers, v, marker=marker, color=color, ms=4, lw=1, label=k)
if ref is not None:
    ax.plot(np.ravel(np.column_stack([ref_lo, ref_hi])), np.repeat(np.where(np.isfinite(ref_v) & (ref_v > 0), ref_v, np.nan), 2),
            color="grey", ls="--", lw=1.5, label="paper stereo 30°")
ax.set(xscale="log", xlim=(0.1, 1000), ylim=(0, 0.5), xlabel=r"$E_\mathrm{true}$ (TeV)", ylabel=r"$\theta_{68}$ (deg)")
ax.grid(True, alpha=0.3)
ax.legend(frameon=False, fontsize=8, title=f"{label}, no cuts", title_fontsize=8, ncol=2)
fig.savefig(out / f"stereo_direction_{label}.png", dpi=150, bbox_inches="tight")
