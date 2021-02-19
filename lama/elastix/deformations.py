#!/usr/bin/env python


"""
Given a sequence of deformation fields, generate a mean deformation field and jacobian determinant file
"""

from lama import common
from logzero import logger as logging
import os
import sys
import subprocess
from pathlib import Path
from typing import Union, Dict, List
import shutil
import SimpleITK as sitk
import numpy as np
from lama.registration_pipeline.validate_config import LamaConfig

ELX_TFORM_NAME = 'TransformParameters.0.txt'
ELX_TFORM_NAME_RESOLUTION = 'TransformParameters.0.R{}.txt'  # resoltion number goes in '{}'
TRANSFORMIX_LOG = 'transformix.log'


def make_deformations_at_different_scales(config: Union[LamaConfig, dict]) -> Union[None, np.array]:
    """
    Generate jacobian determinants and optionaly defromation vectors

    Parameters
    ----------
    config:
        LamaConfig object if running from other lama module
        Path to config file if running this module independently

    Notes
    -----
    How to generate the hacobian determinants is defined by the config['generate_deformation_fields'] entry.

    toml representation from the LAMA config:

    [generate_deformation_fields]
    192_to_10 = [ "deformable_192_to_10",]

    This will create a set of jacobian determinants and optional deformation fields called 192_to_10 and using the
    the named regisrtation stages in the list.

    Multiple sets of key/value pairs are allowed so that diffrent jacobians determinatnts might be made. For eaxmple
    you may want to include the affine transformation in the jacobians, which would look like so:

    affine_192_to_10 = [ "affine", "deformable_192_to_10"]

    Returns
    -------
    jacobian array if there are any negative values
    """

    if isinstance(config, (str, Path)):
        config = LamaConfig(Path(config))

    if not config['generate_deformation_fields']:
        return

    deformation_dir = config.mkdir('deformations')
    jacobians_dir = config.mkdir('jacobians')
    log_jacobians_dir = config.mkdir('log_jacobians')

    write_vectors = config['write_deformation_vectors']
    write_raw_jacobians = config ['write_raw_jacobians']
    write_log_jacobians = config['write_log_jacobians']

    for deformation_id, stage_info in config['generate_deformation_fields'].items():
        reg_stage_dirs: List[Path] = []

        resolutions: List[int] = []  # Specify the resolutions from a stage to generate defs and jacs from

        if len(stage_info) > 1:

            try:
                int(stage_info[1])
            except ValueError:
                # 1. We are spcifiying multiple lama-specified stages: Use the TransformParameters.O.txt from each stage
                for stage_id in stage_info:
                    reg_stage_dirs.append(config['root_reg_dir'] / stage_id)

            else:
                # 2. It is one stage but we are using the elastix defined stages, which will be specified by ints
                #   Use TransformParameters.O.R<x>.txt files where x is the numbers specified
                reg_stage_dirs.append(config['root_reg_dir'] / stage_info[0])
                resolutions = stage_info[1:]
        else:
            # Just one stage is defined so use the TransformParameters.O.txt from that stage
            reg_stage_dirs.append(config['root_reg_dir'] / stage_info[0])

        deformation_id = str(deformation_id)

        deformation_scale_dir = deformation_dir / deformation_id
        deformation_scale_dir.mkdir()
        jacobians_scale_dir =  jacobians_dir / deformation_id
        jacobians_scale_dir.mkdir()
        log_jacobians_scale_dir = log_jacobians_dir / deformation_id
        log_jacobians_scale_dir.mkdir()

        neg_jac_arr = _generate_deformation_fields(reg_stage_dirs, resolutions, deformation_scale_dir, jacobians_scale_dir,
                                                   log_jacobians_scale_dir, write_vectors, write_raw_jacobians, write_log_jacobians,
                                                   threads=config['threads'], filetype=config['filetype'])
        return neg_jac_arr


def _generate_deformation_fields(registration_dirs: List,
                                 resolutions: List,
                                 deformation_dir: Path,
                                 jacobian_dir: Path,
                                 log_jacobians_dir: Path,
                                 write_vectors: bool,
                                 write_raw_jacobians: bool,
                                 write_log_jacobians: bool,
                                 threads=None,
                                 filetype='nrrd',
                                 jacmat=False):
    """
    Run transformix on the specified registration stage to generate deformation fields and spatial jacobians
    """
    logging.info('### Generating deformation files ###')

    specimen_list = [x for x in registration_dirs[0].iterdir() if (registration_dirs[0] / x).is_dir()]

    # if len(specimen_list) < 1:
    #     logging.warn('Can't find any )

    for specimen_path in specimen_list:
        specimen_id = specimen_path.name
        temp_transform_files_dir = deformation_dir / specimen_id
        temp_transform_files_dir.mkdir(exist_ok=True)

        transform_params = []
        # Get the transform parameters for the subsequent registrations

        if len(resolutions) == 0:  # Use the whole stages by using the joint transform file

            for reg_dir in registration_dirs:
                single_reg_dir = reg_dir / specimen_id
                elastix_tform_file = single_reg_dir / ELX_TFORM_NAME  # Transform file in the registration directory

                if not os.path.isfile(elastix_tform_file):
                    raise FileNotFoundError(f"Error. Cannot find elastix transform parameter file: {elastix_tform_file}")

                # Create a new path to copy the transform file to. Then it's in the same forlder as the jacobians etc
                temp_transform_file = temp_transform_files_dir / f'{reg_dir.name}_{specimen_id}_{elastix_tform_file.name}'

                shutil.copy(elastix_tform_file, temp_transform_file)
                transform_params.append(temp_transform_file)
            _modfy_tforms(transform_params)  # Add the InitialtransformParamtere line

        else:
            # The resolutdeformation_dirion paramter files are numbered from 0 but the config counts from 1
            for i in resolutions:
                i -= 1
                single_reg_dir = registration_dirs[0] / specimen_id
                elastix_tform_file = single_reg_dir / ELX_TFORM_NAME_RESOLUTION.format(i)

                if not os.path.isfile(elastix_tform_file):
                    raise FileNotFoundError(f"### Error. Cannot find elastix transform parameter file: {elastix_tform_file}")

                temp_transform_file = temp_transform_files_dir / (registration_dirs[0].name + '_' + elastix_tform_file.name)

                shutil.copy(elastix_tform_file, temp_transform_file)
                transform_params.append(temp_transform_file)

            _modfy_tforms(transform_params)  # Add the InitialtransformParamtere line

        # Copy the tp files into the temp directory and then modify to add initail transform

        # pass in the last tp file [-1] as the other tp files are internally referenced withinn this file
        neg_jac_array = _get_deformations(transform_params[-1], deformation_dir, jacobian_dir, log_jacobians_dir, filetype, specimen_id,
                                          threads, jacmat, write_vectors, write_raw_jacobians, write_log_jacobians)
        return neg_jac_array


def _modfy_tforms(tforms: List):
    """
    Add the initial paramter file paths to the tform files
    Wedon't use this now as all the transforms are merged into one by elastix
    :return:
    """
    if len(tforms) < 2:  # Cannot have initial tform file as we need at least 2
        return

    for i, tp in enumerate(tforms[1:]):
        initial_tp = tforms[i]

        with open(tp, 'r') as fh:
            lines = []
            for line in fh:
                if line.startswith('(InitialTransformParametersFileName'):
                    continue
                lines.append(line)

        previous_tp_str = '(InitialTransformParametersFileName "{}")'.format(initial_tp)
        lines.insert(0, previous_tp_str + '\n')
        with open(tp, 'w') as wh:
            for line in lines:
                wh.write(line)


def _get_deformations(tform: Path,
                      deformation_dir: Path,
                      jacobian_dir: Path,
                      log_jacobians_dir: Path,
                      filetype: str,
                      specimen_id: str,
                      threads: int,
                      make_jacmat: bool,
                      write_vectors: bool = False,
                      write_raw_jacobians: bool = False,
                      write_log_jacobians: bool = True) -> Union[None, np.array]:
    """
    Generate spatial jacobians and optionally deformation files.

    Returns
    -------
    the jacobian array if there are any values < 0
    """

    cmd = ['transformix',
           '-out', str(deformation_dir),
           '-tp', str(tform),
           '-jac', 'all'
           ]

    if write_vectors:
        cmd.extend(['-def', 'all'])
    if make_jacmat:
        cmd.extend(['-jacmat', 'all'])
    if threads:
        cmd.extend(['-threads', str(threads)])

    try:
        subprocess.check_output(cmd)

    except subprocess.CalledProcessError as e:
        logging.exception('transformix failed')
        logging.exception(e)
        # raise subprocess.CalledProcessError(f'### Transformix failed ###\nError message: {e}\nelastix command:{cmd}')
        raise ValueError
    else:
        deformation_out = deformation_dir / f'deformationField.{filetype}'
        jacobian_out = deformation_dir / f'spatialJacobian.{filetype}'

        # rename and move output
        if write_vectors:
            new_def = deformation_dir / (specimen_id + '.' + filetype)
            shutil.move(deformation_out, new_def)

        new_jac = jacobian_dir / (specimen_id + '.' + filetype)

        try:
            shutil.move(jacobian_out, new_jac)
        except IOError:
            #  Bit of a hack. If trasforms conatain subtransforms from pairwise, elastix is unable to generate
            # deformation fields. So try with itk
            def_img = sitk.ReadImage(new_def)
            jac_img = sitk.DisplacementFieldJacobianDeterminant(def_img)
            sitk.WriteImage(jac_img, new_jac)

        # if we have full jacobian matrix, rename and remove that
        if make_jacmat:
            make_jacmat.mkdir()
            jacmat_file = deformation_dir / f'fullSpatialJacobian.{filetype}'  # The name given by elastix
            jacmat_new = make_jacmat / (specimen_id + '.' + filetype)           # New informative name
            shutil.move(jacmat_file, jacmat_new)

        # test if there has been any folding in the jacobians
        jac_img = sitk.ReadImage(str(new_jac))
        jac_arr = sitk.GetArrayFromImage(jac_img)
        jac_min = jac_arr.min()
        jac_max = jac_arr.max()
        logging.info("{} spatial jacobian, min:{}, max:{}".format(specimen_id, jac_min, jac_max))

        if jac_min <= 0:
            logging.warning(
                "The jacobian determinant for {} has negative values. You may need to add a penalty term to the later registration stages".format(
                    specimen_id))
            # Highlight the regions folding
            jac_arr[jac_arr > 0] = 0
            log_jac_path = log_jacobians_dir / ('ERROR_NEGATIVE_JACOBIANS_' + specimen_id + '.' + filetype)
            common.write_array(jac_arr, log_jac_path)

        elif write_log_jacobians:
            # Spit out the log transformed jacobians
            log_jac = np.log(jac_arr)
            log_jac_path = log_jacobians_dir / ( 'log_jac_' + specimen_id + '.' + filetype)

            if not write_raw_jacobians:
                new_jac.unlink()

            common.write_array(log_jac, log_jac_path)


    logging.info('Finished generating deformation fields')

    if jac_min <=0:
        return jac_arr



