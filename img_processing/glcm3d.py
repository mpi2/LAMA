#!/usr/bin/python

"""
Implement a 3D GLCM for use in the registration  pipeline

"""

import numpy as np
import SimpleITK as sitk
import multiprocessing
from multiprocessing import Process
import tempfile
from radiomics import firstorder, getTestCase, glcm, glrlm, glszm, imageoperations, shape
from os.path import join, basename, splitext, dirname, realpath
import uuid
import common
import subprocess
import logging
import yaml

MAXINTENSITY = 255
# GLCM constants
CHUNK_SIZE = 5
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


def pyradiomics_glcm(vol_dir, output_dir, mask, chunksize=10, metric='Contrast'):
    settings = {'binWidth': 25,
            'interpolator': sitk.sitkBSpline,
            'resampledPixelSpacing': None}

    vol_paths = common.get_file_paths(vol_dir)

    for path in vol_paths:
        array = common.img_path_to_array(path)
        result = []

        for chunk in common.get_chunks(array, chunksize, mask):

            chunk_img = sitk.GetImageFromArray(chunk)  # Need to make a sitk chunker
            mask = np.ones_like(chunk)  # No mask. Set all to 1s
            mask_img = sitk.GetImageFromArray(mask)

            glcmFeatures = glcm.RadiomicsGLCM(chunk_img, mask_img, **settings)
            glcmFeatures.enableAllFeatures()
            glcmFeatures.calculateFeatures()
            result.append(glcmFeatures.featureValues[metric])

        # Write the 1D array result. Only where chunks where mask does not all == 0.
        # Can be rebuilt using common.rebuid_subsamlped_output
        out_array = np.array(result)
        out_path = join(out_dir, splitext(basename(path))[0] + '.npy')
        np.save(out_path, out_array)


def itk_glcm_generation(vol_dir, output_dir, chunksize=5, feature_type='all'):
    """

    Parameters
    ----------
    vol_dir: str
        Directory containing volumes. Can be be in sub folders
    output_dir: str
        Where to put the output. Folder must exist
    chunksize: int
        the size of the chunck to make each glcm from
    feature_type: str
        what feature type to report

    Returns
    -------

    """
    feature_types = ["Energy", "Entropy", "InverseDifferenceMoment", "Inertia", "ClusterShade", "ClusterProminence"]

    vol_paths = common.get_file_paths(vol_dir)

    for feature in feature_types:

        feature_out_dir = join(output_dir, feature)
        common.mkdir_if_not_exists(feature_out_dir)

        first = True
        for im_path in vol_paths:
            if first:
                first = False
                i = sitk.ReadImage(im_path)
                a = sitk.GetArrayFromImage(i)
                shape = a.shape
            base = splitext(basename(im_path))[0]
            glcm_outpath = join(feature_out_dir, base + '.bin')

            try:
                out = subprocess.check_output([PATH_TO_ITK_GLCM, im_path, glcm_outpath, feature])
                print out + '\n'
            except subprocess.CalledProcessError as e:
                logging.warn("glcm generation failed: {}".format(e))
                return
            except OSError as e:
                logging.warn("Cannot find path to the texture analysis binary. Has it been compiled for your system: {}".format(e))


        out_config = {
            'original_shape': list(shape),
            'chunksize': chunksize}

        out_config_path = join(feature_out_dir, 'glcm.yaml')

        with open(out_config_path, 'w') as gf:
            gf.write(yaml.dump(out_config))


if __name__ == '__main__':
    import sys
    input_ = sys.argv[1]
    out_dir = sys.argv[2]
    mask = sys.argv[3]

    pyradiomics_glcm(input_, out_dir, mask, chunksize=10, metric='Contrast')





