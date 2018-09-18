#! /usr/bin/env python

"""staging_metric_maker.py
After registration, this module gets a proxy of developmental stage for each input.

The following methods are implemented:
    
    scaling factor (sf)
        The affine or Similarity scaling factors of of each input is then used to determine the CRL
        
        How to use: 
        need: CRL of population average target in um
    
    Whole embryo volume (WEV) - Not yet implemeted
        A tight mask of the population average target is first created
        After registration the mask is backpropagated onto te inputs to get a voxel count of them
    
    Label length (LL) - not yet implemented
        This method was used as it can be more robust than the other two methods. First segment an organ from the 
        population average. We have been using the spine. After registration, this label is back propagated onto the
        inputs. The length of organ is then used as a stage proxy.
"""

from . import affine_similarity_scaling_factors as asf
from .skeleton_length import run as skeleton
from os.path import basename, join, splitext
import os
import common
import fnmatch
from logzero import logger as logging
import pandas as pd

HEADER = 'vol,value\n'


def scaling_factor_staging(root_registration_dir, outdir):
    """
    Make a csv of estimated CRL (from registration scaling factor) viven a list of registration folders
    (usually the final non-linear stage)

    Parameters
    ----------
    root_registration_dir: str
        Path to registration folder containing subfolders of registrations from which to extract scaling factor
    outdir:
    Returns
    -------

    """
    output = {}

    registration_folders = os.listdir(root_registration_dir)
    for dir_ in registration_folders:
        dir_ = join(root_registration_dir, dir_)
        if not os.path.isdir(dir_):
            continue
        tform_param = asf.extract_affine_transformation_parameters(dir_)
        if not tform_param:
            print(dir_)
        scaling_factor = asf.get_scaling_factor(tform_param)
        vol_id = basename(dir_)
        output[vol_id] = scaling_factor
    _write_output(output, outdir)


def whole_volume_staging(inverted_mask_dir, outdir):
    """
    Generate a csv of whole embryo volumes.

    Parameters
    ----------
    inverted_mask_dir: str
        path to folder containing masks that have been inverted back to rigid or original inputs
    """

    output = {}

    inv_mask_folders = os.listdir(inverted_mask_dir)
    for mask_folder_name in inv_mask_folders:
        mask_folder_path = join(inverted_mask_dir, mask_folder_name)
        if not os.path.isdir(mask_folder_path):
            continue
        mask_file_name = fnmatch.filter(os.listdir(mask_folder_path), mask_folder_name + '*')[0]
        mask_path = join(mask_folder_path, mask_file_name)

        mask_array = common.LoadImage(mask_path).array
        embryo_vol_voxels = mask_array[mask_array == 1].size
        output[mask_folder_name] = embryo_vol_voxels
    _write_output(output, outdir)


def label_length_staging(label_inversion_dir, outdir):
    lengths = skeleton(label_inversion_dir)
    _write_output(lengths, outdir)


def _write_output(outdict, outdir):
    outfile = join(outdir, common.STAGING_INFO_FILENAME)
    with open(outfile, 'w') as fh:
        fh.write(HEADER)
        for id_, value in outdict.items():
            fh.write("{},{}\n".format(id_, value))


if __name__=='__main__':
    import sys
    inv_mask_root_dir = sys.argv[1]
    outpath = sys.argv[2]
    whole_volume_staging(inv_mask_root_dir, outpath)