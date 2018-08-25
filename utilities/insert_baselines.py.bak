#! /usr/bin/env python


"""
insert_baselines.py

The following scenario was in mind when making this script.

You have a set of baseline data that has been registered towards your current population average.
It contains, amongst other things:

    registered images
    jabobians
    and possibly inverted labels - which might be needed for calculating staging or organ volumes

This module will insert the above data into the correct locations of the WT test set hierarchy
If you have inverted labels that are used to calculate staging metrics, these will need to be run manually after insertion

example source csv
----
/home/neil/bit/LAMA_results/E14.5/paper_runs/mutant_runs/analysed_lines/sort1/,20150417_SORT1_E14.5_13.3g_WT_XX_REC_scaled_4.6878_pixel_14.0
/home/neil/bit/LAMA_results/E14.5/paper_runs/mutant_runs/analysed_lines/TOR1AIP1/,20170510_TOR1AIP1_E14.5_1.3d_WT_XY_REC_scaled_4.6823_pixel_14.0001
/home/neil/bit/LAMA_results/E14.5/paper_runs/mutant_runs/analysed_lines/TOR1AIP1/,20170511_TOR1AIP1_E14.5_3.1d_WT_XY_rec_scaled_4.6823_pixel_14.0001
/home/neil/bit/LAMA_results/E14.5/paper_runs/mutant_runs/analysed_lines/SLC17A6/,20161124_SLC17A6_E14.5_11.2a_WT_XY_REC_scaled_4.7297_pixel_13.9999
/home/neil/bit/LAMA_results/E14.5/paper_runs/mutant_runs/analysed_lines/SLC17A6/,20160905_SLC17A6_E14.5_8.1d_WT_XX_Rec_scaled_4.7297_pixel_13.999
----
"""

import paths
from os.path import split, join, abspath, isfile, isdir, basename, splitext
import os
import yaml
import common
import shutil

def insert(src_csv, wt_test_set_config_path):
    """

    Parameters
    ----------
    new_config_path: str
        path to config for new baselines
    wt_test_set_config_path: str
        path to config for original baselines
    """

    EXT = '.nrrd'

    with open(src_csv, 'r') as fh:
        for line in fh:
            root, id_ = line.strip().split(',')

            testset_config = yaml.load(open(wt_test_set_config_path, 'r'))
            testset_config_dir = split(abspath(wt_test_set_config_path))[0]

            new_reg_dir = join(root, 'output', 'registrations')
            testset_reg_dir = join(testset_config_dir, 'output', 'registrations')


            ################copy the registrations
            registration_stages = [x['stage_id'] for x in testset_config.get('registration_stage_params')]
            # Get the new reg folders

            for stage in registration_stages:
                new_stage_root = join(new_reg_dir, stage)

                source = join(new_stage_root, id_)
                destination = join(testset_reg_dir, stage, id_)

                if isdir(destination):
                    print("not copying {}\n folder already exists".format(destination))
                else:
                    shutil.copytree(source, destination)
                    print('copied {} to {}'.format(source, destination))

            ################ copy the jacobians
            testset_jac_dir = join(testset_config_dir, 'output', 'jacobians')
            new_jac_dir = join(root, 'output', 'jacobians')
            jacobian_scales = os.listdir(new_jac_dir)  # a scale being 192-16, or 128-16 for example
            for scale in jacobian_scales:
                source_jac_scale_dir = join(new_jac_dir, scale)
                jac_src = join(source_jac_scale_dir, id_ + '.nrrd')
                jac_dest = join(testset_jac_dir, scale, id_ + '.nrrd')

                if isfile(jac_dest):
                    print("not copying {}\n file already exists")
                else:
                    shutil.copy(jac_src, jac_dest)
                    print('copied {} to {}'.format(jac_src, jac_dest))

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("Insert baselines in WT test set")
    parser.add_argument('-s', '--source_csv', dest='source_csv', help='list of root ', required=True)
    parser.add_argument('-d', '--dest_config', dest='dest_config', help='wt test set lama config file', required=True)
    # parser.add_argument('-i', '--vol_ids', dest='vol_ids', help='csv of volume ids to move',
    #                     required=True)
    args = parser.parse_args()
    insert(args.source_csv, args.dest_config)

    