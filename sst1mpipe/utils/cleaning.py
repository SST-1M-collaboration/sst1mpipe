import logging
from copy import deepcopy

import numpy as np
from ctapipe.containers import (
    CameraHillasParametersContainer,
    CameraTimingParametersContainer,
    ImageParametersContainer,
)
from ctapipe.image import (
    ImageCleaner,
    ImageProcessor,
    apply_time_delta_cleaning,
    number_of_islands,
    tailcuts_clean,
)


def get_only_main_island_mask(geom, cleaning_mask):

    num_islands, island_labels = number_of_islands(geom, cleaning_mask)
    n_pixels_on_island = np.bincount(island_labels)
    # First bin is N pixels in no island, i.e. background
    n_pixels_on_island[0] = 0
    # biggest island
    max_island_label = np.argmax(n_pixels_on_island)
    if max_island_label > 0:
        new_image_mask = island_labels == max_island_label
    else:
        new_image_mask = cleaning_mask
    return new_image_mask


class ImageCleanerSST(ImageCleaner):
    """
    Reclean images using the LST-like cleaning, i.e. dynamic picture/boundary technique + time delta. See
    `ctapipe.image.tailcuts_clean`
    `ctapipe.image.apply_time_delta_cleaning`
    """

    average_charge = 0
    stdev_charge = 0
    nsb_level = 0
    config = 0
    frac_rised = 0
    raised = np.zeros(1296)

    def __call__(
        self, tel_id: int, image: np.ndarray, arrival_times=None
    ) -> np.ndarray:
        """
        Apply SST-1M dynamic cleaning. See `ImageCleaner.__call__()`
        """

        defaults = self.config['telescope_defaults']['tel_0'+str(tel_id)]

        for setting in defaults:
            min_nsb_level = setting['min_nsb_level']
            stdev_scaling = setting['stdev_scaling']
            picture_threshold_pe = setting['picture_threshold_pe']
            boundary_threshold_pe = setting['boundary_threshold_pe']
            min_picture_neighbors = setting['min_picture_neighbors']
            keep_isolated = setting['keep_isolated']
            only_main_island = setting['only_main_island']
            min_time_neighbors = setting['min_time_neighbors']
            time_limit_ns = setting['time_limit_ns']
            if self.nsb_level >= min_nsb_level:
                break

        pic_thr=np.maximum(picture_threshold_pe, self.average_charge + stdev_scaling*self.stdev_charge)
        geom = self.subarray.tel[tel_id].camera.geometry
        try:
            self.frac_rised = sum(pic_thr > picture_threshold_pe)/float(len(pic_thr))
            self.raised += pic_thr > picture_threshold_pe
        except Exception:
            self.frac_rised = 0.0

        mask_tailcuts = tailcuts_clean(
            geom,
            image,
            picture_thresh = pic_thr,
            boundary_thresh = boundary_threshold_pe,
            min_number_picture_neighbors = min_picture_neighbors,
            keep_isolated_pixels = keep_isolated,
        )

        time_delta_cleaning_mask = apply_time_delta_cleaning(
            geom,
            mask = mask_tailcuts,
            arrival_times = arrival_times,
            min_number_neighbors = min_time_neighbors,
            time_limit = time_limit_ns
        )

        if only_main_island:
            return get_only_main_island_mask(geom, time_delta_cleaning_mask)
        else:
            return time_delta_cleaning_mask

    # for testing, might disappear in the future
    def dump(self):
        for ii,r in enumerate(self.raised):
            logging.info("pixel %i raised %i times",ii,r)


class ImageCleanerSST_MC(ImageCleaner):
    """
    Reclean images using the LST-like cleaning, i.e. dynamic picture/boundary technique + time delta. See
    `ctapipe.image.tailcuts_clean`
    `ctapipe.image.apply_time_delta_cleaning`
    """

    nsb_level = 0
    config = 0

    def __call__(
        self, tel_id: int, image: np.ndarray, arrival_times=None
    ) -> np.ndarray:
        """
        Apply SST-1M dynamic cleaning. See `ImageCleaner.__call__()`
        """

        defaults = self.config['telescope_defaults']['tel_00'+str(tel_id)]

        for setting in defaults:
            min_nsb_level = setting['min_nsb_level']
            picture_threshold_pe = setting['picture_threshold_pe']
            boundary_threshold_pe = setting['boundary_threshold_pe']
            min_picture_neighbors = setting['min_picture_neighbors']
            keep_isolated = setting['keep_isolated']
            only_main_island = setting['only_main_island']
            min_time_neighbors = setting['min_time_neighbors']
            time_limit_ns = setting['time_limit_ns']
            if self.nsb_level >= min_nsb_level:
                break

        geom = self.subarray.tel[tel_id].camera.geometry

        mask_tailcuts = tailcuts_clean(
            geom,
            image,
            picture_thresh = picture_threshold_pe,
            boundary_thresh = boundary_threshold_pe,
            min_number_picture_neighbors = min_picture_neighbors,
            keep_isolated_pixels = keep_isolated,
        )

        time_delta_cleaning_mask = apply_time_delta_cleaning(
            geom,
            mask = mask_tailcuts,
            arrival_times = arrival_times,
            min_number_neighbors = min_time_neighbors,
            time_limit = time_limit_ns
        )

        if only_main_island:
            return get_only_main_island_mask(geom, time_delta_cleaning_mask)
        else:
            return time_delta_cleaning_mask


def image_cleaner_setup(subarray=None, config=None, ismc=False):

    cleaner = config['ImageProcessor']['image_cleaner_type']

    # If we use other than a standard cleaner recognized by ctapipe
    # (['TailcutsImageCleaner', 'MARSImageCleaner', 'FACTImageCleaner', 'TimeConstrainedImageCleaner'])
    # we have to setup its configuration manualy. Is there a better way to tackle this?
    if cleaner == 'ImageCleanerSST':
        image_processor   = ImageProcessor(subarray=subarray)
        image_processor.use_telescope_frame = config['ImageProcessor']['use_telescope_frame']
        if ~image_processor.use_telescope_frame:
            DEFAULT_IMAGE_PARAMETERS_CAMFRAME = deepcopy(ImageParametersContainer())
            DEFAULT_IMAGE_PARAMETERS_CAMFRAME.hillas = CameraHillasParametersContainer()
            DEFAULT_IMAGE_PARAMETERS_CAMFRAME.timing = CameraTimingParametersContainer()
            image_processor.default_image_container = DEFAULT_IMAGE_PARAMETERS_CAMFRAME

        if ismc:
            image_cleaner_sst = ImageCleaner.from_name(
                    cleaner+'_MC', subarray=subarray, parent=image_processor
                )
            image_processor.clean = image_cleaner_sst
        else:
            image_cleaner_sst = ImageCleaner.from_name(
                    cleaner, subarray=subarray, parent=image_processor
                )
            image_processor.clean = image_cleaner_sst
            image_cleaner_sst.average_charge = 0
            image_cleaner_sst.stdev_charge = 0
            image_cleaner_sst.nsb_level = 0

        image_cleaner_sst.config = config['ImageProcessor'][cleaner]

    else:
        image_processor   = ImageProcessor(subarray=subarray, config=config)

    return image_processor
