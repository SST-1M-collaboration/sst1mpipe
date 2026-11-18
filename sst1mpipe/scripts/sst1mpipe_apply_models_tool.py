"""
``sst1mpipe-apply-models``: ``ctapipe-apply-models`` with the models of the node of a production of
models closest to the pointing of the input file. The output is the DL2 of ctapipe; the models can use
the ``StereoMARSCombiner``.

    sst1mpipe-apply-models -i run.dl1.h5 --models-directory models -o run.dl2.h5
"""
import re
from pathlib import Path as PathlibPath

import astropy.units as u
import numpy as np
from ctapipe.core.traits import Float, Path
from ctapipe.io import TableLoader
from ctapipe.tools.apply_models import ApplyModels as CtapipeApplyModels

import sst1mpipe.ml  # noqa: F401 (StereoMARSCombiner)

__all__ = ["ApplyModels", "find_model_nodes", "closest_model_node", "input_pointing"]

#: subdirectories of a production of models: zenith, azimuth (deg), NSB rate (MHz)
MODEL_NODE_PATTERN = re.compile(r"^ze(?P<zenith>[-\d.]+)_az(?P<azimuth>[-\d.]+)_nsb(?P<nsb>[-\d.]+)$")


def find_model_nodes(models_directory):
    """``{subdirectory: (zenith, azimuth, nsb)}`` of the subdirectories ``zeXX_azXX_nsbXX``"""
    nodes = {}
    for path in sorted(PathlibPath(models_directory).iterdir()):
        match = MODEL_NODE_PATTERN.match(path.name)
        if path.is_dir() and match:
            nodes[path] = tuple(float(match[key]) for key in ("zenith", "azimuth", "nsb"))
    return nodes


def closest_model_node(nodes, zenith, azimuth=None, nsb_rate=None):
    """The node closest in NSB rate (if given), then in zenith, then in azimuth (if given)"""
    if not nodes:
        raise ValueError("No model node (subdirectory zeXX_azXX_nsbXX)")

    def distance(node):
        ze, az, nsb = node[1]
        return (
            0 if nsb_rate is None else abs(nsb - nsb_rate),
            abs(ze - zenith),
            0 if azimuth is None else abs((az - azimuth + 180) % 360 - 180),
        )

    return min(nodes.items(), key=distance)[0]


def input_pointing(path):
    """Median zenith and circular mean azimuth (deg) of the pointing of the telescope events of a DL1 file"""
    with TableLoader(path) as loader:
        events = loader.read_telescope_events(dl1_parameters=False, pointing=True)
    if "telescope_pointing_altitude" in events.colnames:
        altitude, azimuth = events["telescope_pointing_altitude"], events["telescope_pointing_azimuth"]
    else:
        altitude, azimuth = events["subarray_pointing_lat"], events["subarray_pointing_lon"]
    altitude = u.Quantity(altitude, u.deg).to_value(u.rad)
    azimuth = u.Quantity(azimuth, u.deg).to_value(u.rad)
    mean_azimuth = np.arctan2(np.nanmean(np.sin(azimuth)), np.nanmean(np.cos(azimuth)))
    return 90 - np.rad2deg(np.nanmedian(altitude)), np.rad2deg(mean_azimuth) % 360


class ApplyModels(CtapipeApplyModels):
    """
    ``ctapipe-apply-models`` with the models (``*.pkl``) of the node of ``models_directory`` (one
    subdirectory per node ``zeXX_azXX_nsbXX``: zenith, azimuth in deg, NSB rate in MHz) closest to
    the pointing of the input file: closest NSB rate (if ``nsb_rate`` is given), then zenith, then
    azimuth. With ``--reconstructor``, the models given are used instead.
    """

    name = "sst1mpipe-apply-models"

    models_directory = Path(
        default_value=None,
        allow_none=True,
        exists=True,
        file_ok=False,
        help="Production of models: one subdirectory zeXX_azXX_nsbXX per node, with its models (*.pkl)",
    ).tag(config=True)

    nsb_rate = Float(
        default_value=None,
        allow_none=True,
        help="NSB rate of the input file (MHz) to choose the models, not used if not given",
    ).tag(config=True)

    aliases = {
        **CtapipeApplyModels.aliases,
        ("m", "models-directory"): "ApplyModels.models_directory",
        "nsb-rate": "ApplyModels.nsb_rate",
    }

    def setup(self):
        if not self.reconstructor_paths and self.models_directory is not None and self.input_url is not None:
            zenith, azimuth = input_pointing(self.input_url)
            node = closest_model_node(find_model_nodes(self.models_directory), zenith, azimuth, self.nsb_rate)
            self.reconstructor_paths = sorted(node.glob("*.pkl"))
            self.log.info(
                "Pointing zenith %.1f deg, azimuth %.1f deg: models of %s %s",
                zenith, azimuth, node.name, [p.name for p in self.reconstructor_paths],
            )
        super().setup()


def main():
    ApplyModels().run()
