#!/usr/bin/env python
"""
Pointing of the SST-1M telescopes from the logs of their drive system, while tracking on source:

- right ascension and declination of each telescope as a function of time, one column per field
  of the night (successive tracking periods closer than ``--field-radius``)
- map of the time spent tracking in each direction of the sky (2D histogram): whole sky, and
  zoom on each field for each telescope

The pointing given by the TARGET field of zfits runs (``--zfits-files``) is drawn for comparison.

Usage:

$> sst1mpipe_plot_drive_pointing
--input-files DriveSystem1_20260120_001.fits DriveSystem2_20260120_002.fits
--zfits-files SST1M1_20260120_1179.fits.fz
--output-dir ./plots
"""
import argparse
import logging
from pathlib import Path

import astropy.units as u
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord, UnitSphericalRepresentation, concatenate
from astropy.io import fits
from astropy.time import Time
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from astropy.wcs import WCS
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter

from sst1mpipe.io.drive import drive_radec, read_drive_logs
from sst1mpipe.io.sst1m_event_source import parse_target_field

# telescopes: categorical slots 1 and 2 of the reference palette (blue, orange)
TEL_COLORS = {21: "#2a78d6", 22: "#eb6834"}
TEXT = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"
# time spent per bin: one hue, light to dark (blue ramp, steps 200 to 700)
TIME_CMAP = LinearSegmentedColormap.from_list(
    "tracking_time",
    ["#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"],
)
TIME_CMAP.set_bad(alpha=0)

# longest time between two samples of the drive counted as tracking time (s)
MAX_SAMPLE_DURATION = 5


def parse_args():
    parser = argparse.ArgumentParser(description="Pointing of the telescopes from the drive logs")
    parser.add_argument("--input-files", "-f", nargs="+", required=True, help="Drive logs (DriveSystem*.fits)")
    parser.add_argument(
        "--zfits-files", "-z", nargs="+", default=[], help="zfits runs whose TARGET pointing is drawn for comparison"
    )
    parser.add_argument("--output-dir", "-o", default=".", help="Directory of the plots")
    parser.add_argument(
        "--field-radius", type=float, default=10.0,
        help="Successive tracking periods closer than this (deg) are in the same field",
    )
    parser.add_argument(
        "--min-tracking", type=float, default=30.0, help="Tracking periods shorter than this (s) are not used"
    )
    return parser.parse_args()


def mean_direction(coords):
    """Mean direction of sky coordinates"""
    mean = coords.cartesian.mean().represent_as(UnitSphericalRepresentation)
    return SkyCoord(mean.lon, mean.lat, frame=coords.frame.replicate_without_data())


def wrap_ra(ra, center):
    """Right ascension (deg) within 180 deg of center (deg), continuous through RA 0"""
    return center + (ra - center + 180) % 360 - 180


def tracking_periods(drive, min_tracking):
    """
    Samples of the drive logs tracking on source (without the settling at the start of the
    tracking) with their right ascension and declination, split in tracking periods (list of
    index arrays of the samples) at least min_tracking long (s)
    """
    index = np.flatnonzero(drive["is_tracking"] & drive["is_on_source"])
    track = drive[index]
    radec = drive_radec(track)
    # without the AltAz attributes (time, location) carried by the SkyCoord, to merge the telescopes
    track["radec"] = SkyCoord(radec.ra, radec.dec, frame=radec.frame.replicate_without_data())
    duration = np.diff(track["time"].unix, append=np.inf)
    track["duration"] = np.clip(duration, 0, MAX_SAMPLE_DURATION) * u.s
    periods = np.split(np.arange(len(track)), np.flatnonzero(np.diff(index) > 1) + 1)
    periods = [p for p in periods if (track["time"][p[-1]] - track["time"][p[0]]).to_value(u.s) >= min_tracking]
    return track, periods


def find_fields(tracks, radius):
    """
    Fields of the night: tracking periods of all the telescopes, in time order, closer than
    radius to the first period of the field. Each field has its time range (start, stop) and
    its center (mean direction of its periods).
    """
    periods = []
    for track, track_periods in tracks.values():
        for p in track_periods:
            periods.append((track["time"][p[0]], track["time"][p[-1]], mean_direction(track["radec"][p])))
    periods.sort(key=lambda period: period[0])

    fields = []
    for start, stop, center in periods:
        if fields and fields[-1]["first"].separation(center) < radius:
            fields[-1]["stop"] = max(fields[-1]["stop"], stop)
            fields[-1]["centers"].append(center)
        else:
            fields.append(dict(start=start, stop=stop, first=center, centers=[center]))
    for field in fields:
        field["center"] = mean_direction(concatenate(field["centers"]))
    return fields


def read_targets(files):
    """Name, start time and pointing of the TARGET field of zfits runs (runs without pointing are skipped)"""
    targets = []
    for file in files:
        with fits.open(file) as hdul:
            header = hdul["Events"].header
        target, wobble, ra, dec = parse_target_field(header.get("TARGET"))
        if ra is None:
            logging.warning("No pointing in the TARGET field (%r) of %s", header.get("TARGET"), file)
            continue
        targets.append(dict(
            target=target,
            name=target if wobble in (None, "UNDEF") else f"{target} {wobble}",
            time=Time(header["DATE"], scale="utc"),
            coord=SkyCoord(ra * u.deg, dec * u.deg, frame="icrs"),
        ))
    return targets


def field_targets(field, targets, radius):
    """Targets of the zfits runs in the field, one per pointing"""
    selected = {}
    for target in targets:
        if field["start"] <= target["time"] <= field["stop"] and target["coord"].separation(field["center"]) < radius:
            selected[target["name"]] = target
    return list(selected.values())


def field_name(field, targets):
    names = sorted({target["target"] for target in targets})
    if names:
        return ", ".join(names)
    center = field["center"]
    return f"RA {center.ra.deg:.1f}°, Dec {center.dec.deg:+.1f}°"


def in_field(track, field):
    return (track["time"] >= field["start"]) & (track["time"] <= field["stop"])


def style_axes(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=8, length=3)


def plot_radec_vs_time(tracks, fields, targets, radius, path):
    """Right ascension and declination of the telescopes while tracking, one column per field"""
    widths = [max((f["stop"] - f["start"]).to_value(u.h), 0.5) for f in fields]
    fig, axes = plt.subplots(
        2, len(fields), figsize=(15, 6.5), sharex="col", squeeze=False,
        gridspec_kw=dict(width_ratios=widths, hspace=0.08, wspace=0.22),
    )
    fig.patch.set_facecolor(SURFACE)
    ra_formatter = FuncFormatter(lambda value, _: f"{value % 360:g}")

    for column, field in enumerate(fields):
        ax_ra, ax_dec = axes[:, column]
        center_ra = field["center"].ra.deg
        for tel_id, (track, periods) in tracks.items():
            mask = in_field(track, field)
            for p in periods:
                p = p[mask[p]]
                if len(p) == 0:
                    continue
                time = track["time"][p].plot_date
                radec = track["radec"][p]
                kwargs = dict(color=TEL_COLORS.get(tel_id, TEXT_SECONDARY), linewidth=1.4, solid_capstyle="round")
                ax_ra.plot(time, wrap_ra(radec.ra.deg, center_ra), **kwargs)
                ax_dec.plot(time, radec.dec.deg, **kwargs)

        for target in field_targets(field, targets, radius):
            coord = target["coord"].transform_to(field["center"].frame)
            for ax, value in ((ax_ra, wrap_ra(coord.ra.deg, center_ra)), (ax_dec, coord.dec.deg)):
                ax.axhline(value, color=TEXT_SECONDARY, linewidth=1, linestyle=(0, (4, 3)))
                ax.annotate(
                    f"{target['name']} (TARGET)", (1, value), xycoords=("axes fraction", "data"),
                    xytext=(-2, 3), textcoords="offset points", ha="right", va="bottom",
                    fontsize=7.5, color=TEXT_SECONDARY,
                )

        ax_ra.set_title(field_name(field, field_targets(field, targets, radius)), fontsize=10, color=TEXT, loc="left")
        ax_ra.yaxis.set_major_formatter(ra_formatter)
        for ax in (ax_ra, ax_dec):
            style_axes(ax)
            ax.margins(y=0.12)
        ax_dec.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=2, maxticks=6))
        ax_dec.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

    axes[0, 0].set_ylabel("RA (deg)", color=TEXT_SECONDARY, fontsize=9)
    axes[1, 0].set_ylabel("Dec (deg)", color=TEXT_SECONDARY, fontsize=9)
    fig.supxlabel("Time (UTC)", color=TEXT_SECONDARY, fontsize=9, y=0.02)

    night = min(f["start"] for f in fields).to_datetime().strftime("%Y-%m-%d")
    fig.suptitle(f"Pointing of the drive while tracking on source, night of {night}", x=0.06, ha="left", fontsize=13, color=TEXT)
    handles = [plt.Line2D([], [], color=TEL_COLORS.get(tel_id, TEXT_SECONDARY), linewidth=2) for tel_id in tracks]
    labels = [f"Tel {tel_id}" for tel_id in tracks]
    if targets:
        handles.append(plt.Line2D([], [], color=TEXT_SECONDARY, linewidth=1, linestyle=(0, (4, 3))))
        labels.append("TARGET of the zfits run")
    fig.legend(handles, labels, loc="upper right", ncols=len(handles), frameon=False, fontsize=9, labelcolor=TEXT)
    fig.savefig(path, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    logging.info("Plot saved in %s", path)


def allsky_wcs(fields, frame):
    """
    Hammer-Aitoff projection of the whole sky (1 deg pixels, the sky is an ellipse of half axes
    2 sqrt(2) and sqrt(2) rad), with the fields away from its edge
    """
    ras = np.sort([field["center"].ra.deg for field in fields])
    gaps = np.diff(np.append(ras, ras[0] + 360))
    # the edge of the map (center + 180 deg) in the middle of the largest gap between the fields
    edge = ras[np.argmax(gaps)] + gaps.max() / 2
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---AIT", "DEC--AIT"]
    wcs.wcs.crval = [(edge + 180) % 360, 0]
    wcs.wcs.cdelt = [-1, 1]
    half_x, half_y = int(np.ceil(np.degrees(2 * np.sqrt(2)))), int(np.ceil(np.degrees(np.sqrt(2))))
    wcs.wcs.crpix = [half_x + 1, half_y + 1]
    wcs.wcs.radesys = frame.name.upper()
    wcs.wcs.equinox = frame.equinox.jyear
    return wcs, (2 * half_y + 1, 2 * half_x + 1)


def field_wcs(center, half_width, n_bins=100):
    """Gnomonic projection centered on the field, n_bins x n_bins pixels"""
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [center.ra.deg, center.dec.deg]
    scale = 2 * half_width / n_bins
    wcs.wcs.cdelt = [-scale, scale]
    wcs.wcs.crpix = [(n_bins + 1) / 2, (n_bins + 1) / 2]
    wcs.wcs.radesys = center.frame.name.upper()
    wcs.wcs.equinox = center.frame.equinox.jyear
    return wcs, (n_bins, n_bins)


def histogram(wcs, shape, radec, weights):
    """Sum of the weights in the pixels of the map (NaN in the empty pixels)"""
    x, y = wcs.world_to_pixel(radec)
    counts, _, _ = np.histogram2d(
        y, x, bins=(np.arange(shape[0] + 1) - 0.5, np.arange(shape[1] + 1) - 0.5), weights=weights
    )
    counts[counts == 0] = np.nan
    return counts


def galactic_plane(wcs, frame):
    """Pixel coordinates of the Galactic plane, NaN where it crosses the edge of the map"""
    plane = SkyCoord(np.linspace(0, 360, 721) * u.deg, np.zeros(721) * u.deg, frame="galactic").transform_to(frame)
    x, y = wcs.world_to_pixel(plane)
    jumps = np.flatnonzero(np.hypot(np.diff(x), np.diff(y)) > 20) + 1
    return np.insert(x, jumps, np.nan), np.insert(y, jumps, np.nan)


def style_map(ax, ra_spacing=None, dec_spacing=None, label_ra=True, label_dec=True):
    ax.set_facecolor(SURFACE)
    ra, dec = ax.coords[0], ax.coords[1]
    ax.coords.grid(color=GRID, linewidth=0.6, linestyle="solid")
    for coord, spacing, label, name in ((ra, ra_spacing, label_ra, "RA (deg)"), (dec, dec_spacing, label_dec, "Dec (deg)")):
        coord.set_major_formatter("d.d" if spacing is None or spacing < 1 * u.deg else "d")
        if spacing is not None:
            coord.set_ticks(spacing=spacing, color=GRID)
        coord.set_ticklabel(color=TEXT_SECONDARY, size=7.5)
        coord.set_axislabel(name if label else "", color=TEXT_SECONDARY, size=8.5)
        coord.set_ticklabel_visible(label)
    ax.coords.frame.set_color(GRID)


def nice_spacing(half_width):
    """Spacing of the grid lines of a map of half width half_width (deg)"""
    for spacing in (0.25, 0.5, 1, 2, 5):
        if 2 * half_width / spacing <= 6:
            return spacing * u.deg
    return 10 * u.deg


def plot_radec_map(tracks, fields, targets, radius, path):
    """Time spent tracking in each direction: whole sky and zoom on each field for each telescope"""
    frame = fields[0]["center"].frame.replicate_without_data()
    tel_ids = list(tracks)
    n_columns = len(fields)
    fig = plt.figure(figsize=(max(4.6 * n_columns, 10), 3.6 + 4.2 * len(tel_ids) + 0.6))
    fig.patch.set_facecolor(SURFACE)
    grid = fig.add_gridspec(
        2 + len(tel_ids), n_columns, height_ratios=[3.6] + [4.2] * len(tel_ids) + [0.1], hspace=0.25, wspace=0.12,
        top=0.95, bottom=0.04,
    )

    # whole sky
    wcs, shape = allsky_wcs(fields, frame)
    ax = fig.add_subplot(grid[0, :], projection=wcs, frame_class=EllipticalFrame)
    radec = concatenate([track["radec"] for track, _ in tracks.values()])
    minutes = np.concatenate([track["duration"].to_value(u.min) for track, _ in tracks.values()])
    image = ax.imshow(histogram(wcs, shape, radec, minutes), cmap=TIME_CMAP, origin="lower", interpolation="nearest")
    ax.plot(*galactic_plane(wcs, frame), color=TEXT_SECONDARY, linewidth=0.8, alpha=0.6)
    x, y = wcs.world_to_pixel(SkyCoord(100 * u.deg, 0 * u.deg, frame="galactic").transform_to(frame))
    ax.text(x + 3, y, "Galactic plane", fontsize=7, color=TEXT_SECONDARY, va="center")
    style_map(ax, ra_spacing=60 * u.deg, dec_spacing=30 * u.deg)
    ax.coords[0].set_ticklabel_visible(False)
    ax.coords[0].set_axislabel("")
    ax.coords[1].set_axislabel("")
    for value in range(0, 360, 60):
        x, y = wcs.world_to_pixel(SkyCoord(value * u.deg, 0 * u.deg, frame=frame))
        ax.text(x, y - 3, f"{value}°", ha="center", va="top", fontsize=7, color=TEXT_SECONDARY)
    for field in fields:
        x, y = wcs.world_to_pixel(field["center"])
        ax.annotate(
            field_name(field, field_targets(field, targets, radius)), (x, y), xytext=(8, 8), textcoords="offset points",
            fontsize=8.5, color=TEXT, arrowprops=dict(arrowstyle="-", color=TEXT_SECONDARY, linewidth=0.6),
        )
    colorbar = fig.colorbar(image, cax=ax.inset_axes([1.0, 0.15, 0.012, 0.7]))
    colorbar.set_label("Tracking time per 1° pixel (min)", color=TEXT_SECONDARY, fontsize=8)
    colorbar.ax.tick_params(colors=TEXT_SECONDARY, labelsize=7.5)
    colorbar.outline.set_visible(False)
    ax.set_title("All telescopes, whole sky (Hammer-Aitoff, FK5)", fontsize=10, color=TEXT, loc="left")

    # zoom on each field, one row per telescope
    for column, field in enumerate(fields):
        selected = field_targets(field, targets, radius)
        samples = {tel_id: in_field(track, field) for tel_id, (track, _) in tracks.items()}
        separations = np.concatenate([
            tracks[tel_id][0]["radec"][mask].separation(field["center"]).deg for tel_id, mask in samples.items()
        ])
        half_width = max(1.15 * np.max(separations, initial=0), 0.5)
        wcs, shape = field_wcs(field["center"], half_width)
        images = {
            tel_id: histogram(wcs, shape, tracks[tel_id][0]["radec"][mask], tracks[tel_id][0]["duration"][mask].to_value(u.min))
            for tel_id, mask in samples.items()
        }
        vmax = max(np.nanmax(image, initial=0) for image in images.values())
        dec_spacing = nice_spacing(half_width)
        ra_spacing = nice_spacing(half_width / np.cos(field["center"].dec.rad))
        for row, tel_id in enumerate(tel_ids, start=1):
            ax = fig.add_subplot(grid[row, column], projection=wcs)
            image = ax.imshow(images[tel_id], cmap=TIME_CMAP, origin="lower", vmin=0, vmax=vmax, interpolation="nearest")
            for target in selected:
                x, y = wcs.world_to_pixel(target["coord"].transform_to(frame))
                ax.plot(x, y, marker="+", markersize=11, markeredgewidth=1.4, color=TEXT)
                ax.annotate(
                    target["name"], (x, y), xytext=(6, -10), textcoords="offset points", fontsize=7.5, color=TEXT
                )
            # each field has its own coordinates: RA labels on the last row, Dec labels on all the maps
            style_map(ax, ra_spacing, dec_spacing, label_ra=row == len(tel_ids))
            if column > 0:
                ax.coords[1].set_axislabel("")
            ax.text(
                0.03, 0.96, f"Tel {tel_id}", transform=ax.transAxes, ha="left", va="top", fontsize=9, color=TEXT,
                bbox=dict(boxstyle="round,pad=0.25", facecolor=SURFACE, edgecolor=TEL_COLORS.get(tel_id, GRID)),
            )
            if row == 1:
                ax.set_title(field_name(field, selected), fontsize=10, color=TEXT, loc="left")
        cax = fig.add_subplot(grid[-1, column])
        colorbar = fig.colorbar(image, cax=cax, orientation="horizontal")
        bin_size = 2 * half_width / shape[0] * 60
        colorbar.set_label(f"Tracking time per {bin_size:.1f}′ bin (min)", color=TEXT_SECONDARY, fontsize=8)
        colorbar.ax.tick_params(colors=TEXT_SECONDARY, labelsize=7.5)
        colorbar.outline.set_visible(False)

    night = min(f["start"] for f in fields).to_datetime().strftime("%Y-%m-%d")
    fig.suptitle(f"Sky map of the drive pointing while tracking on source, night of {night}", x=0.06, y=0.995, ha="left", fontsize=13, color=TEXT)
    fig.savefig(path, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    logging.info("Plot saved in %s", path)


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tracks = {}
    for file in args.input_files:
        drive = read_drive_logs(file)
        tracks[drive.meta["tel_id"]] = tracking_periods(drive, args.min_tracking)
    tracks = dict(sorted(tracks.items()))
    radius = args.field_radius * u.deg
    fields = find_fields(tracks, radius)
    targets = read_targets(args.zfits_files)
    for field in fields:
        logging.info(
            "Field %s from %s to %s", field_name(field, field_targets(field, targets, radius)),
            field["start"].isot, field["stop"].isot,
        )

    plot_radec_vs_time(tracks, fields, targets, radius, output_dir / "drive_radec_vs_time.png")
    plot_radec_map(tracks, fields, targets, radius, output_dir / "drive_radec_map.png")


if __name__ == "__main__":
    main()
