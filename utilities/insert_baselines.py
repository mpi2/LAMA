"""
insert_baselines.py

The following scenario was in mind when making this scipt.

You have a set of baseline data that has been registered towards your current population average.
It contains, amongst other things:

    registered images
    jabobians
    and possibly inverted labels - which might be needed for calculating staging or organ volumes

This module will insert the above data into the correct locations of the WT test set hierarchy
If you have inverted labels that are used to calculate staging metrics, these will need to be run manually after insertion
"""

import paths
from os.path import split, join, abspath, isfile, isdir, basename, splitext
import os
import yaml
import common
import shutil

def insert(new_config_path, wt_test_set_config_path):
    """

    Parameters
    ----------
    new_config_path: str
        path to config for new baselines
    wt_test_set_config_path: str
        path to config for original baselines
    """

    new_config = yaml.load(open(new_config_path, 'r'))
    new_config_dir = split(abspath(new_config_path))[0]
    new_reg_paths = paths.RegPaths(new_config_dir, new_config)
    new_data_input_dir = new_reg_paths.get('inputs')
    vol_ids = common.strip_img_extensions(common.get_file_paths(new_data_input_dir))

    testset_config = yaml.load(open(wt_test_set_config_path, 'r'))
    testset_config_dir = split(abspath(wt_test_set_config_path))[0]
    testset_paths = paths.RegPaths(testset_config_dir, testset_config)

    ################copy the registrations and inverted labels
    testset_reg_dir = testset_paths.get('registrations')
    new_reg_dir = new_reg_paths.get('registrations')
    registration_stages = [x['stage_id'] for x in new_config.get('registration_stage_params')]
    # Get the new reg folders

    for stage in registration_stages:
        new_stage_root = join(new_reg_dir, stage)

        sources = [join(new_stage_root, x) for x in vol_ids]
        destinations = [join(testset_reg_dir, stage, x) for x in vol_ids]

        for src, dst in zip(sources, destinations):
            if isdir(dst):
                print("not copying {}\n folder already exists".format(dst))
            else:
                shutil.copytree(src, dst)
                print('copied {} to {}'.format(src, dst))

    ################ copy the jacobians
    testset_jac_dir = testset_paths.get('jacobians')
    new_jac_dir = new_reg_paths.get('jacobians')
    jacobian_scales = os.listdir(new_jac_dir)  # a scale being 192-16, or 128-16 for example
    for scale in jacobian_scales:
        source_jac_scale_dir = join(new_jac_dir, scale)
        source_filenames = [basename(x) for x in os.listdir(source_jac_scale_dir) if splitext(basename(x))[0] in vol_ids]
        j_dests = [join(testset_jac_dir, scale, x) for x in source_filenames]
        j_sources = [join(new_jac_dir, scale, x) for x in source_filenames]
        for src, dst in zip(j_sources, j_dests):
            if isfile(dst):
                print("not copying {}\n file already exists")
            else:
                shutil.copy(src, dst)
                print('copied {} to {}'.format(src, dst))

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("Insert baselines in WT test set")
    parser.add_argument('-s', '--source_config', dest='source_config', help='source lama config file', required=True)
    parser.add_argument('-d', '--dest_config', dest='dest_config', help='wt test set lama config file', required=True)
    # parser.add_argument('-i', '--vol_ids', dest='vol_ids', help='csv of volume ids to move',
    #                     required=True)
    args = parser.parse_args()
    insert(args.source_config, args.dest_config)

    