"""
``sst1mpipe-dl2-dl3``: DL3 of the data (GADF, for gammapy) from the DL2 of ``sst1mpipe-apply-models``,
with the cuts and point-like IRFs of the ctapipe IRF tools (``ctapipe-optimize-event-selection``,
``ctapipe-compute-irf``) of the zenith angle closest to the observation.

    sst1mpipe-dl2-dl3 -i stereo/run.dl2.h5 --irf-directory irf/optimized/stereo -o dl3/

One DL3 file per observation (``<output-dir>/SST1M_<obs_id>.dl3.fits.gz``) in the Gamma Astro Data
Formats (GADF, https://gamma-astro-data-formats.readthedocs.io, see `sst1mpipe.io.dl3`), then the index
files of all the DL3 files of the output directory (``hdu-index.fits.gz``, ``obs-index.fits.gz``):

    from gammapy.data import DataStore
    observations = DataStore.from_dir("dl3/").get_observations(required_irf="point-like")
"""
import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, SkyCoord
from astropy.io import fits
from astropy.table import Table, join
from astropy.time import Time
from ctapipe.core import Tool
from ctapipe.core.traits import Bool, Path, Unicode, flag
from ctapipe.io import TableLoader
from ctapipe.irf import OptimizationResult

from sst1mpipe.io.dl3 import (
    events_hdu,
    event_selection,
    gti_hdu,
    irf_hdus,
    pointing_hdu,
    select_irf_directory,
    write_index,
)
from sst1mpipe.scripts.sst1mpipe_apply_models_tool import add_telescope_pointing

__all__ = ["DL2ToDL3"]


def _circular_mean(angles):
    """Mean of angles (deg) on the circle"""
    radians = np.deg2rad(angles)
    return np.rad2deg(np.arctan2(np.mean(np.sin(radians)), np.mean(np.cos(radians)))) % 360


def _target(observation):
    """ICRS position of the target of the observation (observation block), None if not in ICRS"""
    from ctapipe.containers import CoordinateFrameType

    names = ("subarray_pointing_frame", "subarray_pointing_lon", "subarray_pointing_lat")
    if not all(name in observation.colnames for name in names):
        return None
    try:
        frame = CoordinateFrameType(int(observation["subarray_pointing_frame"][0]))
    except (ValueError, TypeError):
        return None
    lon, lat = observation["subarray_pointing_lon"][0], observation["subarray_pointing_lat"][0]
    if frame is not CoordinateFrameType.ICRS or not np.isfinite([float(lon), float(lat)]).all():
        return None
    return SkyCoord(ra=u.Quantity(lon, u.deg), dec=u.Quantity(lat, u.deg), frame="icrs")


class DL2ToDL3(Tool):
    """DL3 of the data from the DL2 of ctapipe, with the cuts and IRFs of the ctapipe IRF tools"""

    name = "sst1mpipe-dl2-dl3"
    description = __doc__
    examples = "sst1mpipe-dl2-dl3 -i run.dl2.h5 --irf-directory irf/optimized/stereo -o dl3/"

    input_url = Path(exists=True, directory_ok=False, help="DL2 file of sst1mpipe-apply-models").tag(config=True)
    output_dir = Path(file_ok=False, help="Directory of the DL3 files and of the index files").tag(config=True)
    irf_directory = Path(
        exists=True, file_ok=False,
        help=(
            "Directory with cuts.fits and irf_point_like.fits.gz (ctapipe-optimize-event-selection,"
            " ctapipe-compute-irf), or with such directories ze<zenith in deg>: the closest zenith angle is used"
        ),
    ).tag(config=True)
    geometry_reconstructor = Unicode("disp", help="Prefix of the reconstructed direction").tag(config=True)
    energy_reconstructor = Unicode("RandomForestRegressor", help="Prefix of the reconstructed energy").tag(config=True)
    telescope = Unicode("SST-1M", help="TELESCOP of the DL3 files").tag(config=True)
    origin = Unicode("SST-1M Collaboration", help="ORIGIN of the DL3 files (organisation)").tag(config=True)
    overwrite = Bool(False, help="Overwrite the DL3 files").tag(config=True)
    index = Bool(
        True, help="Write the index files of all the DL3 files of the output directory (not when run in parallel)",
    ).tag(config=True)

    aliases = {
        ("i", "input"): "DL2ToDL3.input_url",
        ("o", "output-dir"): "DL2ToDL3.output_dir",
        "irf-directory": "DL2ToDL3.irf_directory",
    }
    flags = {
        **flag("overwrite", "DL2ToDL3.overwrite", "Overwrite the DL3 files", "Do not overwrite the DL3 files"),
        **flag("index", "DL2ToDL3.index", "Write the index files", "Do not write the index files"),
    }

    def setup(self):
        if self.input_url is None or self.output_dir is None or self.irf_directory is None:
            raise ValueError("--input, --output-dir and --irf-directory are required")
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def _events(self, loader):
        """Subarray events with their pointing (alt/az in deg)"""
        events = loader.read_subarray_events(dl2=True, observation_info=True)
        pointing = None
        try:
            telescopes = loader.read_telescope_events(dl1_parameters=False, dl2=False, pointing=True,
                                                      observation_info=False)
            if "telescope_pointing_altitude" in telescopes.colnames:
                pointing = telescopes
        except Exception as error:  # noqa: BLE001 (no pointing table)
            self.log.debug("No pointing of the telescopes: %s", error)
        if pointing is not None:
            # mean pointing of the telescopes of each event (azimuth: on the circle)
            az = np.deg2rad(u.Quantity(pointing["telescope_pointing_azimuth"], u.deg).value)
            per_event = Table({
                "obs_id": pointing["obs_id"], "event_id": pointing["event_id"],
                "pointing_alt": u.Quantity(pointing["telescope_pointing_altitude"], u.deg).value,
                "sin_az": np.sin(az), "cos_az": np.cos(az),
            }).group_by(["obs_id", "event_id"]).groups.aggregate(np.mean)
            per_event["pointing_az"] = np.rad2deg(np.arctan2(per_event["sin_az"], per_event["cos_az"])) % 360
            per_event.remove_columns(["sin_az", "cos_az"])
            events = join(events, per_event, keys=["obs_id", "event_id"], join_type="left")
            events.sort(["obs_id", "time"])
            for name in ("pointing_alt", "pointing_az"):
                events[name] = np.ma.filled(np.ma.asarray(events[name], dtype=float), np.nan)
        else:
            self.log.warning("No pointing of the telescopes in %s: pointing of the target used", self.input_url)
            if not add_telescope_pointing(events, loader.subarray.reference_location):
                raise ValueError(f"No pointing in {self.input_url}")
            events["pointing_alt"] = u.Quantity(events["telescope_pointing_altitude"], u.deg).value
            events["pointing_az"] = u.Quantity(events["telescope_pointing_azimuth"], u.deg).value
        return events

    def start(self):
        with TableLoader(self.input_url) as loader:
            subarray = loader.subarray
            events = self._events(loader)
        location = subarray.reference_location
        geometry, energy = self.geometry_reconstructor, self.energy_reconstructor
        self.dl3_files = []
        for obs_id in np.unique(events["obs_id"]):
            observation = events[events["obs_id"] == obs_id]
            zenith = 90 - np.nanmean(observation["pointing_alt"])
            irf_directory = select_irf_directory(self.irf_directory, zenith * u.deg)
            cuts = OptimizationResult.read(irf_directory / "cuts.fits")
            gammaness = f"{cuts.clf_prefix}_prediction"
            selected = observation[event_selection(observation, cuts, f"{energy}_energy", gammaness)]
            self.log.info("obs_id %d: %d / %d events, IRFs of %s (zenith %.1f deg)", obs_id, len(selected),
                          len(observation), irf_directory, zenith)

            time = Time(observation["time"])
            start, stop = time.min(), time.max()
            dl3_events = {
                "event_id": selected["event_id"], "time": selected["time"],
                "energy": u.Quantity(selected[f"{energy}_energy"], u.TeV),
                "alt": u.Quantity(selected[f"{geometry}_alt"], u.deg), "az": u.Quantity(selected[f"{geometry}_az"], u.deg),
                "gammaness": selected[gammaness],
            }
            if f"{geometry}_telescopes" in selected.colnames:
                # telescopes of the reconstruction of the direction (MULTIP)
                dl3_events["multiplicity"] = np.sum(np.asarray(selected[f"{geometry}_telescopes"]), axis=1)
            middle = start + (stop - start) / 2
            alt, az = np.nanmean(observation["pointing_alt"]), _circular_mean(observation["pointing_az"])
            pointing = SkyCoord(alt=alt * u.deg, az=az * u.deg, frame=AltAz(obstime=middle, location=location)).icrs
            import sst1mpipe

            tel_ids = [int(t) for t in subarray.tel_ids]
            meta = {
                "ORIGIN": (self.origin, "Organisation that created the file"),
                "TELESCOP": (self.telescope, "Telescope"),
                # SST-1M-1, SST-1M-2 (data telescopes 21, 22, simulated 1, 2) or SST-1M-1+2 (stereo)
                "INSTRUME": ("SST-1M-" + "+".join(str(t % 20) for t in tel_ids), "Telescopes of the analysis"),
                "N_TELS": (len(tel_ids), "Number of telescopes"),
                "TELLIST": (",".join(str(t) for t in tel_ids), "Telescope ids"),
                "ALT_PNT": (alt, "[deg] Mean pointing altitude"), "AZ_PNT": (az, "[deg] Mean pointing azimuth"),
                "CREATOR": (f"{self.name} (sst1mpipe {sst1mpipe.__version__})", "Software that created the file"),
                "IRF": str(irf_directory),  # IRFs and cuts (ctapipe-compute-irf); long path: no comment
            }
            if "target" in observation.colnames and str(observation["target"][0]):
                meta["OBJECT"] = (str(observation["target"][0]), "Target of the observation")
            target = _target(observation)
            if target is not None:
                meta["RA_OBJ"] = (target.ra.to_value(u.deg), "[deg] Right ascension of the target")
                meta["DEC_OBJ"] = (target.dec.to_value(u.deg), "[deg] Declination of the target")
            hdus = [
                events_hdu(dl3_events, obs_id, location, pointing, start, stop, meta),
                gti_hdu(start, stop),
                pointing_hdu(time, observation["pointing_alt"], observation["pointing_az"], location),
                *irf_hdus(irf_directory / "irf_point_like.fits.gz"),
            ]
            path = self.output_dir / f"SST1M_{obs_id}.dl3.fits.gz"
            fits.HDUList([fits.PrimaryHDU(), *hdus]).writeto(path, overwrite=self.overwrite)
            self.dl3_files.append(path)

    def finish(self):
        if not self.index:
            return
        dl3_files = sorted(self.output_dir.glob("*.dl3.fits.gz"))
        write_index(dl3_files, self.output_dir)
        self.log.info("%d DL3 files indexed in %s", len(dl3_files), self.output_dir)


def main():
    DL2ToDL3().run()
