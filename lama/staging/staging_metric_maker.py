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
from pathlib import Path
from os.path import join
import os
from typing import Dict

from lama.staging import affine_similarity_scaling_factors as asf
from lama.staging.skeleton_length import run as skeleton

from lama import common

from lama import seg_methods
from seg_methods import mask

HEADER = 'vol,value\n'

DEFAULT_STAGING_METHOD = 'embryo_volume'


def scaling_factor_staging(root_registration_dir: Path, outdir: Path):
    """
    Make a csv of estimated CRL (from registration scaling factor) viven a list of registration folders
    (usually the final non-linear stage)

    Parameters
    ----------
    root_registration_dir: Contains subfolders of registrations from which to generate staging data from
    outdir: Whaere to put staging csv
    """

    root_registration_dir = Path(root_registration_dir)
    outdir = Path(outdir)

    output = {}

    for dir_ in root_registration_dir.iterdir():

        if not os.path.isdir(dir_):
            continue

        tform_param = asf.extract_affine_transformation_parameters(dir_)

        if not tform_param:
            print(dir_)

        scaling_factor = asf.get_scaling_factor(tform_param)
        vol_id = dir_.stem
        output[vol_id] = scaling_factor

    _write_output(output, outdir)
    return scaling_factor


def whole_volume_staging(inverted_mask_dir: Path, outdir: Path):
    """
    Generate a csv of whole embryo volumes.

    Parameters
    ----------
    inverted_mask_dir:  masks that have been inverted back to rigid or original inputs
    outdir: where to put the resulting staging csv

    """
    output = {}

    for mask_folder in inverted_mask_dir.iterdir():

        if not mask_folder.is_dir():
            continue

        # The mask with start with the same name as the folder + an image extension
        mask_path = common.getfile_startswith(mask_folder, mask_folder.name)

        mask_array = common.LoadImage(mask_path).array
        embryo_vol_voxels = mask_array[mask_array == 1].size
        output[mask_folder.name] = embryo_vol_voxels
    _write_output(output, outdir)


def label_length_staging(label_inversion_dir, outdir):
    lengths = skeleton(label_inversion_dir)
    _write_output(lengths, outdir)


def otsu_staging(original_vol_dir: Path, outdir: Path, otsudir_root: Path):
    """
    generates a csv file of WEVs labelled by the otsu method

    Parameters
    ----------
    Orginal vol dir
    outdir
        Dir to output staging info

    """
    outdir = Path(outdir)
    otsudir_root = Path(otsudir_root)
    original_vol_dir = Path(original_vol_dir)

    output = {}
    otsudir_root.mkdir(exist_ok=True)

    # Compute otsu for each image. If running job_runner there will only be one image per folder
    # But if running lama_reg there may be multiple

    for dir_ in original_vol_dir.iterdir():

        if not dir_.is_dir():
            continue

        input_image_file = dir_ / f"{dir_.name}.nrrd"

        otsu_specimen_dir = otsudir_root / dir_.name
        otsu_specimen_dir.mkdir(exist_ok=True)
        arguments = ['-tight','-external','-option=otsu','-ndilate','-i{}'.format(input_image_file),'-o{}'.format(otsu_specimen_dir)]
        mask.main(arguments) #generate the otsu volume  # Not how you should call python code
        
        #The results from otsu are a bit weird and need renaming,
        #this gets the mask we want (otsu_binary_img.nrrd), moves it up a directory and renames it something useful.

        # Neil: maybe edit the otsu code to give proprely named outputs. It's not my code so can't help there

        print(os.listdir(otsu_specimen_dir))
        
        #otsufile = Path(otsudir) / "otsu_binary_img.nrrd"
        
        #for result in os.listdir(otsudir): 
        #    if "otsu_binary" in result: 
        #        result = Path(result).absolute()
        #        parent_dir = result.parents[1]
        #        result.rename(parent_dir / "whole_mask_{}".format(file))
        #    else: 
        #        os.remove(result)
               
        #Now we need to create the csv files
        mask_name = otsu_specimen_dir / f"{otsu_specimen_dir.name}.nrrd"
        mask_array = common.LoadImage(mask_name).array
        embryo_vol_voxels = mask_array[mask_array == 1].size
        output[otsu_specimen_dir.name] = embryo_vol_voxels

    _write_output(output, outdir)


def _write_output(data: Dict, outdir: Path):
    outfile = outdir / common.STAGING_INFO_FILENAME

    with open(outfile, 'w') as fh:
        fh.write(HEADER)
        for id_, value in data.items():
            fh.write("{},{}\n".format(id_, value))


STAGING_METHODS = {
        'none': lambda x: print('no staging'),
        'otsu': otsu_staging,
        'label_len': label_length_staging,
        'embryo_volume': whole_volume_staging
        # 'scaling_factor': scaling_factor_staging
}

if __name__=='__main__':
    import sys
    inv_mask_root_dir = sys.argv[1]
    outpath = sys.argv[2]
    whole_volume_staging(inv_mask_root_dir, outpath)