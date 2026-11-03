import astropy.units as u
import numpy as np
from ctapipe.core.traits import FloatTelescopeParameter, IntTelescopeParameter
from ctapipe.image import ImageCleaner
from scipy.spatial.distance import cdist
from scipy.sparse import issparse

from sklearn.neighbors import radius_neighbors_graph


def clean_dbscan_fast(sparse_connectivity_matrix, weights, min_points):

    if not issparse(sparse_connectivity_matrix):
        RuntimeWarning(
            "Connectivity matrix is not sparse it might reduce computation performances"
        )

    n_points = sparse_connectivity_matrix.dot(weights)
    mask_core = n_points >= min_points
    mask_border = (sparse_connectivity_matrix.dot(mask_core)) > 0
    mask = mask_core | mask_border

    return mask


class DBSCANImageCleaner(ImageCleaner):
    """
    An image cleaner based on the sklearn.cluster.DBSCAN algorithm
    """

    minimum_pe = IntTelescopeParameter(
        default_value=30, help="Minimum number of p.e. in cluster"
    ).tag(config=True) # TODO make this density per pixel

    picture_threshold_pe = FloatTelescopeParameter(
        default_value=0.0,
        help="Minimum number of p.e. in the image the pixel.",
    ).tag(config=True)

    epsilon_r = FloatTelescopeParameter(
        default_value=38.0, help="Scale parameter for spatial coordinates (in mm)"
    ).tag(config=True)

    def __init__(self, subarray, config=None, parent=None, **kwargs):
        super().__init__(subarray, config, parent, **kwargs)

        self._precompute_distances()

    def __call__(
        self, tel_id: int, image: np.ndarray, arrival_times: np.ndarray = None
    ) -> np.ndarray:

        mask = clean_dbscan_fast(
            self._distances[tel_id],
            weights=image,
            min_points=self.minimum_pe.tel[tel_id],
        )
        mask = mask & (
            image > self.picture_threshold_pe.tel[tel_id]
        )  # negative value are not accepted for Hillas computation

        if (
            mask.sum() <= 1
        ):  # Cleaning with one pixel fails the timing computation (impossible with 0)
            mask[...] = False

        return mask

    def _precompute_distances(self):

        self._distances = {}
        for tel_id in self.subarray.tel_ids:
            geometry = self.subarray.tel[tel_id].camera.geometry
            epsilon_r = self.epsilon_r.tel[tel_id] * u.mm
            x = np.column_stack([geometry.pix_x, geometry.pix_y]) / epsilon_r

            d = radius_neighbors_graph(
                x.to(u.dimensionless_unscaled),
                radius=1.0,
                mode="connectivity",
                include_self=True,
                n_jobs=-1,
            )

            self._distances[tel_id] = d


class TimeDBSCANImageCleaner(ImageCleaner):
    """
    An image cleaner based on the sklearn.cluster.DBSCAN algorithm that uses the image and peak time image
    """

    minimum_pe = IntTelescopeParameter(
        default_value=30, help="Minimum number of p.e. in cluster"
    ).tag(config=True) # TODO make this density per pixel

    picture_threshold_pe = FloatTelescopeParameter(
        default_value=0.0,
        help="Minimum number of p.e. in the image the pixel.",
    ).tag(config=True)

    epsilon_r = FloatTelescopeParameter(
        default_value=38.0, help="Scale parameter for spatial coordinates (in mm)"
    ).tag(config=True)

    epsilon_t = FloatTelescopeParameter(
        default_value=40.0, help="Scale parameter for time (in ns)"
    ).tag(config=True)

    def __init__(self, subarray, config=None, parent=None, **kwargs):
        super().__init__(subarray, config, parent, **kwargs)

        self._precompute_distances_squared()

    def __call__(
        self, tel_id: int, image: np.ndarray, arrival_times: np.ndarray
    ) -> np.ndarray:

        times = arrival_times / self.epsilon_t.tel[tel_id]
        d = (times[:, None] - times[None, :]) ** 2 + self._distances_squared[tel_id]
        d = np.sqrt(d) <= 1.0

        mask = clean_dbscan_fast(
            d, weights=image, min_points=self.minimum_pe.tel[tel_id]
        )
        mask = mask & (image > self.picture_threshold_pe.tel[tel_id])

        if (
            mask.sum() <= 1
        ):  # Cleaning with one pixel fails the timing computation (impossible with 0)
            mask[...] = False

        return mask

    def _precompute_distances_squared(self):
        self._distances_squared = {}

        for tel_id in self.subarray.tel_ids:
            geometry = self.subarray.tel[tel_id].camera.geometry
            x = np.column_stack([geometry.pix_x, geometry.pix_y])
            d = cdist(x.value, x.value) * x.unit
            d /= self.epsilon_r.tel[tel_id] * u.mm

            self._distances_squared[tel_id] = d**2




class DBSCANTimeImageCleaner(ImageCleaner):
    """
    An image cleaner based on the sklearn.cluster.DBSCAN algorithm that uses the peak time image distances weighted
    by the image intensity
    """

    minimum_pe = IntTelescopeParameter(
        default_value=30, help="Minimum number of p.e. in cluster"
    ).tag(config=True) # TODO make this density per pixel

    picture_threshold_pe = FloatTelescopeParameter(
        default_value=0.0,
        help="Minimum number of p.e. in the image the pixel.",
    ).tag(config=True)

    epsilon_t = FloatTelescopeParameter(
        default_value=40.0, help="Scale parameter for time (in ns)"
    ).tag(config=True)

    def __init__(self, subarray, config=None, parent=None, **kwargs):

        super().__init__(subarray, config, parent, **kwargs)

    def __call__(self, tel_id: int, image: np.ndarray, arrival_times: np.ndarray) -> np.ndarray:

        times = arrival_times / self.epsilon_t.tel[tel_id]
        d = np.abs(times[:, None] - times[None, :])
        d = d <= 1.0

        mask = clean_dbscan_fast(
            d, weights=image, min_points=self.minimum_pe.tel[tel_id]
        )
        mask = mask & (image > self.picture_threshold_pe.tel[tel_id])

        if (
            mask.sum() <= 1
        ):  # Cleaning with one pixel fails the timing computation (impossible with 0)
            mask[...] = False

        return mask


class DBSCANImageCleaner3D(ImageCleaner):
    """
    An image cleaner based on the sklearn.cluster.DBSCAN algorithm that uses the waveforms
    """

    minimum_pe = IntTelescopeParameter(
        default_value=30, help="Minimum number of p.e. in cluster"
    ).tag(config=True) # TODO make this density per pixel

    epsilon_r = FloatTelescopeParameter(
        default_value=38.0, help="Scale parameter for spatial coordinates (in mm)"
    ).tag(config=True)

    epsilon_t = FloatTelescopeParameter(
        default_value=40.0, help="Scale parameter for time (in ns)"
    ).tag(config=True)

    min_samples = IntTelescopeParameter(
        default_value=1, help="Minimum number of time samples required"
    ).tag(config=True)

    def __init__(self, subarray, config=None, parent=None, **kwargs):
        super().__init__(subarray, config, parent, **kwargs)

        self._precompute_distances()

    def __call__(
        self, tel_id: int, waveform: np.ndarray, arrival_times=None
    ) -> np.ndarray:

        mask = clean_dbscan_fast(
            self._distances[tel_id],
            weights=waveform.ravel(),
            min_points=self.minimum_pe.tel[tel_id],
        )

        mask = mask.reshape(
            (
                self.subarray.tel[tel_id].camera.readout.n_pixels,
                self.subarray.tel[tel_id].camera.readout.n_samples,
            )
        )
        mask = mask.sum(axis=-1)
        mask = mask >= self.min_samples.tel[tel_id]

        if (
            mask.sum() <= 1
        ):  # Cleaning with one pixel fails the timing computation (impossible with 0)
            return np.zeros(
                self.subarray.tel[tel_id].camera.readout.n_pixels, dtype=bool
            )

        return mask

    def _precompute_distances(self):

        self._distances = {}

        for tel_id in self.subarray.tel_ids:
            geometry = self.subarray.tel[tel_id].camera.geometry
            readout = self.subarray.tel[tel_id].camera.readout

            t = np.arange(readout.n_samples) / readout.sampling_rate

            indices_xyt = np.column_stack(
                [
                    np.repeat(
                        geometry.pix_x.to(u.mm).value / self.epsilon_r.tel[tel_id],
                        readout.n_samples,
                    ),
                    np.repeat(
                        geometry.pix_y.to(u.mm).value / self.epsilon_r.tel[tel_id],
                        readout.n_samples,
                    ),
                    np.tile(
                        t.to(u.ns).value / self.epsilon_t.tel[tel_id], geometry.n_pixels
                    ),
                ]
            )

            d = radius_neighbors_graph(
                indices_xyt, 1.0, mode="connectivity", include_self=True
            )

            self._distances[tel_id] = d
