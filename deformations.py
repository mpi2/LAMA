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
import yaml


ELX_TFORM_NAME = 'TransformParameters.0.txt'
ELX_TFORM_NAME_RESOLUTION = 'TransformParameters.0.R{}.txt'  # resoltion number goes in '{}'
TRANSFORMIX_LOG = 'transformix.log'


def make_deformations_at_different_scales(config_path, root_reg_dir, outdir, threads=4):
    if not isinstance(config_path, dict):
        try:
            config = yaml.load(open(config_path[0], 'r'))
        except Exception as e:
            logging.error("deformations.py cannot open yaml file")
            sys.exit("can't read the YAML config file - {}".format(config_path))
    else:
        config = config_path

    deformation_dir = join(outdir, 'deformations')
    jacobians_dir = join(outdir, 'jacobians')
    log_jacobians_dir = join(outdir, 'log_jacobians')
    common.mkdir_if_not_exists(deformation_dir)
    common.mkdir_if_not_exists(jacobians_dir)
    common.mkdir_if_not_exists(log_jacobians_dir)

    # Example config entry for generating deformation fields
    # This generates deformation/jacobians from the def_128_8 stage using resolutions 1-3
    #
    #generate_deformation_fields:
    #    def_128_8: [affine]
    #    def_64_8: [affine, 2, 3]


    for deformation_id, stage_info in config['generate_deformation_fields'].iteritems():


        if len(stage_info) > 1:  # Not using all the stage - the resolutions to generate defs from have been specified
            resolutions = stage_info[1:]
        else:
            resolutions = []
        stage_id = stage_info[0]

        deformation_id = str(deformation_id)  #
        reg_stage_dir = join(root_reg_dir, stage_id)
        deformation_scale_dir = join(deformation_dir, deformation_id)
        jacobians_scale_dir = join(jacobians_dir, deformation_id)
        log_jacobians_scale_dir = join(log_jacobians_dir, deformation_id)
        common.mkdir_if_not_exists(deformation_scale_dir)
        common.mkdir_if_not_exists(jacobians_scale_dir)
        common.mkdir_if_not_exists(log_jacobians_scale_dir)

        generate_deformation_fields(reg_stage_dir, resolutions, deformation_scale_dir, jacobians_scale_dir,
                                    log_jacobians_scale_dir, threads=threads)


def generate_deformation_fields(registration_dir, resolutions, deformation_dir, jacobian_dir, log_jacobians_dir, threads=None,
                                filetype='nrrd', jacmat=False):
    """
    Run transformix on the specified registration stage to generate deformation fields and spatial jacobians

    Parameters
    ----------
    registration_dirs: list
        registration directories in order that they were registered
    resolutions: list
        list of resolutions from a stage to generate deformations from
    """
    logging.info('### Generating deformation files ###')
    logging.info('Generating deformation fields from following reg_stage {}'.format(registration_dir))

    specimen_list = [x for x in os.listdir(registration_dir) if os.path.isdir(join(registration_dir, x))]

    # if len(specimen_list) < 1:
    #     logging.warn('Can't find any )

    for specimen_id in specimen_list:
        temp_tp_files = join(deformation_dir, specimen_id)
        common.mkdir_if_not_exists(temp_tp_files)

        transform_params = []
        # Get the transform parameters for the subsequent registrations

        if len(resolutions) == 0:  # Use the whole stage by using the joint transform file
            single_reg_dir = join(registration_dir, specimen_id)
            elx_tform_file = join(single_reg_dir, ELX_TFORM_NAME)
            if not os.path.isfile(elx_tform_file):
                sys.exit("### Error. Cannot find elastix transform parameter file: {}".format(elx_tform_file))

            temp_tp_file = join(temp_tp_files, basename(registration_dir) + '_' + basename(elx_tform_file))
            shutil.copy(elx_tform_file, temp_tp_file)
            transform_params.append(temp_tp_file)

        else:
            # The resolution paramter files are numbered from 0 but the config counts from 1
            for i in resolutions:
                i -= 1
                single_reg_dir = join(registration_dir, specimen_id)
                elx_tform_file = join(single_reg_dir, ELX_TFORM_NAME_RESOLUTION.format(i))
                if not os.path.isfile(elx_tform_file):
                    sys.exit("### Error. Cannot find elastix transform parameter file: {}".format(elx_tform_file))

                temp_tp_file = join(temp_tp_files, basename(registration_dir) + '_' + basename(elx_tform_file))
                shutil.copy(elx_tform_file, temp_tp_file)
                transform_params.append(temp_tp_file)
                modfy_tforms(transform_params)  # Add the InitialtransformParamtere line

        # Copy the tp files into the temp directory and then modify to add initail transform


        # pass in the last tp file [-1] as the other tp files are intyernally referenced
        get_deformations(transform_params[-1], deformation_dir, jacobian_dir, log_jacobians_dir, filetype, specimen_id,
                         threads, jacmat)

        # Move the transformix log
        logfile = join(deformation_dir, TRANSFORMIX_LOG)
        shutil.move(logfile, temp_tp_files)

def modfy_tforms(tforms):
    """
    Add the initial paramter file paths to the tform files
    Wedon't use this now as all the transforms are merged into one by elastix
    :return:
    """
    for i, tp in enumerate(tforms[1:]):
        initial_tp = tforms[i]
        with open(tp, 'rb') as fh:
            lines = []
            for line in fh:
                if line.startswith('(InitialTransformParametersFileName'):
                    continue
                lines.append(line)
        previous_tp_str = '(InitialTransformParametersFileName "{}")'.format(initial_tp)
        lines.insert(0, previous_tp_str + '\n')
        with open(tp, 'wb') as wh:
            for line in lines:
                wh.write(line)


def get_deformations(tform, deformation_dir, jacobian_dir, log_jacobians_dir, filetype, specimen_id, threads,
                     jacmat_dir):
    """
    """
    common.mkdir_if_not_exists(log_jacobians_dir) # Delete

    cmd = ['transformix',
           '-out', deformation_dir,
           '-def', 'all',
           '-jac', 'all',
           '-tp', tform
           ]
    if jacmat_dir:
        cmd.extend(['-jacmat', 'all'])
    if threads:
        cmd.extend(['-threads', threads])

    try:
        output = subprocess.check_output(cmd)
    except Exception as e:  # Can't seem to log CalledProcessError
        logging.warn('transformix failed {}'.format(', '.join(cmd)))
        sys.exit('### Transformix failed ###\nError message: {}\nelastix command:{}'.format(e, cmd))
    else:
        deformation_out = join(deformation_dir, 'deformationField.{}'.format(filetype))
        jacobian_out = join(deformation_dir, 'spatialJacobian.{}'.format(filetype))

        # rename and move output
        new_def = join(deformation_dir, specimen_id + '.' + filetype)
        shutil.move(deformation_out, new_def)

        new_jac = join(jacobian_dir, specimen_id + '.' + filetype)
        shutil.move(jacobian_out, new_jac)

        # if we have full jacobian matrix, rename and remove that
        if jacmat_dir:
            common.mkdir_if_not_exists(jacmat_dir)
            jacmat_out = join(deformation_dir, 'fullSpatialJacobian.{}'.format(filetype))
            jacmat_new = join(jacmat_dir, specimen_id + '.' + filetype)
            shutil.move(jacmat_out, jacmat_new)

        #shutil.rmtree(temp_def_dir)

        # test if there has been any folding in the jacobians
        jac_img = sitk.ReadImage(new_jac)
        jac_arr = sitk.GetArrayFromImage(jac_img)
        jac_min = jac_arr.min()
        jac_max = jac_arr.max()
        logging.info("{} spatial jacobian, min:{}, max:{}".format(specimen_id, jac_min, jac_max))
        if jac_min <= 0:
            logging.warn("The jacobian determinant for {} has negative values. You may need to add a penalty term to the later registration stages")
            # Highlight the regions folding
            jac_arr[jac_arr > 0] = 0
            neg_jac = sitk.GetImageFromArray(jac_arr)
            log_jac_path = join(log_jacobians_dir, 'ERROR_NEGATIVE_JACOBIANS_' + specimen_id + '.' + filetype)
            sitk.WriteImage(neg_jac, log_jac_path, True)
        else:
            # Spit out the log transformed jacobians
            log_jac = np.log10(jac_arr)
            log_jac_img = sitk.GetImageFromArray(log_jac)
            log_jac_path = join(log_jacobians_dir, 'log_jac_' + specimen_id + '.' + filetype)
            sitk.WriteImage(log_jac_img, log_jac_path, True)

    logging.info('Finished generating deformation fields')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser("Generate deformation fields and spatial jacobians from elastix registration")
    parser.add_argument('-c', '--config', dest='config', help='LAMA config file', required=True, nargs='*')
    parser.add_argument('-o', '--out', dest='out_dir', help='folder to put results in', required=True)
    parser.add_argument('-r', '--reg_dir', dest='reg_dir', help='directory containing registration output', required=True)
    parser.add_argument('-t', '--threads', dest='threads', help='Numberof threads to use', required=True, type=str)

    args = parser.parse_args()

    make_deformations_at_different_scales(args.config, args.reg_dir, args.out_dir, args.threads)