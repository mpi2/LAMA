#! /usr/bin/ev/python

"""staging.py
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

import affine_scaling_factors as asf
from skeleton_length import run as skeleton
from os.path import basename, join
import os
import common

HEADER = 'vol,value\n'


def scaling_factor_staging(root_registration_dir, outdir):
    """
    Given a list of registration folders (usually the final non-linear stage) get a list of CRLs for each specimen
    Parameters
    ----------
    root_registration_dir
        a list of folders in which to look for elastix registration files
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
            print dir_
        scaling_factor = asf.get_scaling_factor([tform_param])
        vol_id = basename(dir_)
        output[vol_id] = scaling_factor
    outfile = join(outdir, common.STAGING_INFO_FILENAME)
    with open(outfile, 'w') as fh:
        fh.write(HEADER)
        for id_, value in output.iteritems():
            fh.write("{},{}\n".format(id_, value))

def whole_volume_staging():
    pass


def label_length_staging(label_inversion_dir, outdir):
    lengths = skeleton(label_inversion_dir)
    out_path = join(outdir, common.STAGING_INFO_FILENAME)
    with open(out_path, 'w') as fh:
        fh.write(HEADER)
        for id_, length in lengths.iteritems():
            fh.write("{},{}\n".format(id_, length))



if __name__=='__main__':
    path = '/home/neil/work_hdd/work/test_dataset_for_lama_dev/wt/output/inverted_labels/rigid'
    label_length_staging(path, path)