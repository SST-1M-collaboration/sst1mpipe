import numpy as np

from ctapipe.visualization import CameraDisplay

from sst1mpipe.instrument.camera import Camera


GEOMETRY = Camera().geometry
def test_add_colorbar():

    image = np.ones(GEOMETRY.n_pixels)
    display = CameraDisplay(geometry=GEOMETRY, image=image)
    display.add_colorbar()
