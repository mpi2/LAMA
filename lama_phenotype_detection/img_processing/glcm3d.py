#!/usr/bin/python

"""
Implement a 3D GLCM for use in the registration  pipeline

"""

import numpy as np
import SimpleITK as sitk
from logzero import logger as logging
try:
    from radiomics import firstorder, getTestCase, glcm, glrlm, glszm, imageoperations, shape
except ImportError:
    logging.warn('pyradiomics not installed. No glcms will be made')
    pyrad_installed = False
else:
    pyrad_installed = True
from os.path import join, basename, splitext, dirname, realpath
import common
import yaml

MAXINTENSITY = 255
# GLCM constants
CHUNK_SIZE = 10
GLCM_BINS = 8

SCRIPT_DIR = dirname(realpath(__file__))
PATH_TO_ITK_GLCM = join(SCRIPT_DIR, '../dev/texture/GLCMItk/LamaITKTexture')


def _reshape_data(shape, chunk_size, result_data):
    """
    The data from the GLCM analysis is subsampled and so smaller than the original data. To be able to overlay
    onto real image data, we need to upsample the result

    Parameters
    ----------

    Returns
    -------
    A numpy ndarray? Should be 1D
    """
    out_array = np.zeros(shape)
    i = 0
    # We go x-y-z as thats how it comes out of the GLCM generator
    for x in range(0, shape[2] - chunk_size, chunk_size):
        for y in range(0, shape[1] - chunk_size, chunk_size):
            for z in range(0, shape[0] - chunk_size, chunk_size):
                out_array[z: z + chunk_size, y: y + chunk_size, x: x + chunk_size] = result_data[i]
                i += 1

    return out_array


def pyradiomics_glcm(vol_dir, out_dir, mask, chunksize=CHUNK_SIZE, feature='Contrast'):
    """
    Create glcm and xtract features. Spit out features per chunk as a 1D numpy array
    This 1d array can be reassembled into a 3D volume using common.rebuid_subsamlped_output

    Parameters
    ----------
    vol_dir: str
        Directory containing volumes. Can be be in sub folders
    output_dir: str
        Where to put the output. Folder must exist
    mask: numpy.ndarray
        mask used to exclude the masked regions from analysis
    chunksize: int
        the size of the chunck to make each glcm from
    feature: str
        what feature type to report

    """
    if not pyrad_installed:
        return

    settings = {'binWidth': 4,
            'interpolator': sitk.sitkBSpline,
            'resampledPixelSpacing': None}

    vol_paths = common.get_file_paths(vol_dir)

    glcm_mask = np.ones([chunksize] * 3)  # No glcm mask. Set all to 1s. Needed for pyradiomics
    glcm_mask_img = sitk.GetImageFromArray(glcm_mask)

    for path in vol_paths:
        array = common.img_path_to_array(path)
        result = []

        for chunk in common.get_chunks(array, chunksize, mask):

            chunk_img = sitk.GetImageFromArray(chunk)  # I need to make a sitk chunker
            glcmFeatures = glcm.RadiomicsGLCM(chunk_img, glcm_mask_img, **settings)
            glcmFeatures.enableAllFeatures()
            glcmFeatures.calculateFeatures()
            result.append(glcmFeatures.featureValues[feature])

        # Write the 1D array result. Only where chunks where mask does not all == 0.
        out_array = np.array(result)
        out_path = join(out_dir, splitext(basename(path))[0] + '.npy')
        np.save(out_path, out_array)

        out_config = {
            'original_shape': list(array.shape),
            'chunksize': chunksize}

        out_config_path = join(out_dir, 'glcm.yaml')
        with open(out_config_path, 'w') as fh:
            fh.write(yaml.dump(out_config))


if __name__ == '__main__':
    import sys
    input_ = sys.argv[1]
    out_dir = sys.argv[2]
    mask_path = sys.argv[3]
    mask_array = common.img_path_to_array(mask_path)

    pyradiomics_glcm(input_, out_dir, mask_array, feature='Contrast')





