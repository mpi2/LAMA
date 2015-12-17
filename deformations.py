#!/usr/bin/env python


"""
Given a sequence of deformation fields, generate a mean deformation field and jacobian determinant file
"""

import common
import logging
import os
import sys
from os.path import join, basename
import subprocess
import shutil
import SimpleITK as sitk
import numpy as np


def generate_deformation_fields(registration_dirs, deformation_dir, jacobian_dir, filetype='nrrd'):
    """
    Run transformix on the specified registration stage to generate deformation fields and spatial jacobians

    Parameters
    ----------
    registration_dirs: list
        registration directories in order that they were registered
    """
    #logging.info('### Generating deformation files ###')
    # Get the transform parameters

    first_reg_dir = registration_dirs[0]

    specimen_list = [x for x in os.listdir(first_reg_dir) if os.path.isdir(join(first_reg_dir, x))]

    temp_def_dir = join(deformation_dir, 'temp_deformation')
    common.mkdir_if_not_exists(temp_def_dir)
    for specimen_id in specimen_list:
        # Create a tempfile fr storing deformations
        temp_def_files = []

        for i, reg_dir in enumerate(registration_dirs):

            single_reg_dir = join(reg_dir, specimen_id)

            elx_tform_file = join(single_reg_dir, 'TransformParameters.0.txt')
            if not os.path.isfile(elx_tform_file):
                sys.exit("### Error. Cannot find elastix transform parameter file: {}".format(elx_tform_file))

            cmd = ['transformix',
                   '-tp', elx_tform_file,
                   '-out', temp_def_dir,
                   '-def', 'all',
                   ]

            try:
                output = subprocess.check_output(cmd)
            except Exception as e:  # Can't seem to log CalledProcessError
                logging.warn('transformix failed {}'.format(', '.join(cmd)))
                sys.exit('### Transformix failed ###\nError message: {}\nelastix command:{}'.format(e, cmd))
            else:
                # read in and sum up the deformation fields
                deformation_out = join(temp_def_dir, 'deformationField.{}'.format(filetype))

                def_field = sitk.GetArrayFromImage(sitk.ReadImage(deformation_out))
                if i == 0:  #The first def
                    summed_def = def_field
                else:
                    summed_def += def_field
                os.remove(deformation_out)

        mean_def = summed_def / (i + 1) # TODO work on summeddef directly don't create a copy
        mean_def_image = sitk.GetImageFromArray(mean_def)
        sitk.WriteImage(mean_def_image, join(deformation_dir, specimen_id + '.' + filetype), True)

        mean_jac = sitk.DisplacementFieldJacobianDeterminant(mean_def_image)
        sitk.WriteImage(mean_jac, join(jacobian_dir, specimen_id + '.' + filetype), True)

    shutil.rmtree(temp_def_dir)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser("Generate deformation fields and spatial jacobians from elastix registration")
    parser.add_argument('-r', '--reg_dirs', dest='reg_dirs', help='Series of registration directories', required=True, nargs='*')
    parser.add_argument('-d', '--def_out', dest='def_out', help='folder to put deformations in', required=True)
    parser.add_argument('-j', '--jac_out', dest='jac_out', help='folder to put jacobians in', required=True)
    args = parser.parse_args()

    generate_deformation_fields([os.path.abspath(x) for x in args.reg_dirs],
                                os.path.abspath(args.def_out),
                                os.path.abspath(args.jac_out))
