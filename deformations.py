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
import tempfile
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

    for specimen_id in specimen_list:
        # Create a tempfile fr storing deformations
        temp_def_files = []
        temp_def_dir = tempfile.mkdtemp()

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
                # Rename the outputs and move to correct dirs
                deformation_out = join(temp_def_dir, 'deformationField.{}'.format(filetype))
                deformation_renamed = join(temp_def_dir, '{}.{}'.format(str(i), filetype))
                shutil.move(deformation_out, deformation_renamed)
                temp_def_files.append(deformation_renamed)

        # Create total deformations and jacobians
        for c, def_file in enumerate(temp_def_files):
            if c == 0:
                def_array = sitk.GetArrayFromImage(sitk.ReadImage(def_file))
                #shutil.rmtree(def_file)
            else:
                new_def = sitk.GetArrayFromImage(sitk.ReadImage(def_file))
                def_array += new_def
                #shutil.rmtree(def_file)
        mean_def = def_array / (c + 1)
        mean_def_image = sitk.GetImageFromArray(mean_def)
        sitk.WriteImage(mean_def_image, join(deformation_dir, specimen_id + '_deformations.' + filetype), True)

        mean_jac = sitk.DisplacementFieldJacobianDeterminant(mean_def_image)
        sitk.WriteImage(mean_jac, join(jacobian_dir, specimen_id + 'jacobians.' + filetype), True)

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
