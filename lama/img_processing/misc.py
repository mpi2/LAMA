import numpy as np
import SimpleITK as sitk
import numpy as np
from scipy import ndimage



def blur(img: np.ndarray, fwhm: float, voxel_size: float) -> np.ndarray:
    """
    https://matthew-brett.github.io/teaching/random_fields.html

    Parameters
    ----------

    """
    fwhm_in_voxels = fwhm / voxel_size

    sd = fwhm_in_voxels / np.sqrt(8. * np.log(2))  # sigma for this FWHM
    blurred = ndimage.filters.gaussian_filter(img, sd, mode='constant', cval=0.0)

    return blurred
