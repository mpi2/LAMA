from pathlib import Path
import tempfile
import os
import subprocess
from collections import defaultdict
from multiprocessing import Pool
from os.path import join, abspath, isfile
from typing import Union, List, Dict

from logzero import logger as logging
import yaml

from lama import common
from lama.common import cfg_load
from lama.registration_pipeline.validate_config import LamaConfig

from lama.elastix import (ELX_TRANSFORM_NAME, ELX_PARAM_PREFIX, LABEL_INVERTED_TRANFORM,
                          IMAGE_INVERTED_TRANSFORM, INVERT_CONFIG, RESOLUTION_IMGS_DIR, IMG_PYRAMID_DIR)

LABEL_REPLACEMENTS = {
    'FinalBSplineInterpolationOrder': '0',
    'FixedInternalImagePixelType': 'short',
    'MovingInternalImagePixelType': 'short',
    'ResultImagePixelType': 'unsigned char',
    'WriteTransformParametersEachResolution': 'false',
    'WriteResultImageAfterEachResolution': 'false'
}

IMAGE_REPLACEMENTS = {
    'FinalBSplineInterpolationOrder': '3',
    'FixedInternalImagePixelType': 'float',
    'MovingInternalImagePixelType': 'float',
    'ResultImagePixelType': 'float',
    'WriteTransformParametersEachResolution': 'false',
    'WriteResultImageAfterEachResolution': 'false'

}


def batch_invert_transform_parameters(config: Union[str, LamaConfig],
                                      clobber=True, new_log:bool=False):
    """
    Create new elastix TransformParameter files that can then be used by transformix to invert labelmaps, stats etc

    Parameters
    ----------
    config
        path to original reg pipeline config file

    clobber
        if True overwrite inverted parameters present

    new_log:
        Whether to create a new log file. If called from another module, logging may happen there
    """
    common.test_installation('elastix')

    if isinstance(config, (Path, str)):
        config = LamaConfig(config)

    threads = str(config['threads'])

    if new_log:
        common.init_logging(config / 'invert_transforms.log')

    reg_dirs = get_reg_dirs(config)

    # Get the image basenames from the first stage registration folder (usually rigid)
    # ignore images in non-relevent folder that may be present
    volume_names = [x.stem for x in common.get_file_paths(reg_dirs[0], ignore_folders=[RESOLUTION_IMGS_DIR, IMG_PYRAMID_DIR])]

    inv_outdir = config.mkdir('inverted_transforms')

    stages_to_invert = defaultdict(list)

    jobs: List[Dict] = []

    reg_stage_dir: Path

    for i, vol_id in enumerate(volume_names):

        for reg_stage_dir in reg_dirs:

            if not reg_stage_dir.is_dir():
                logging.error('cannot find {}'.format(reg_stage_dir))
                raise FileNotFoundError(f'Cannot find registration dir {reg_stage_dir}')

            inv_stage_dir = inv_outdir / reg_stage_dir.name

            specimen_stage_reg_dir = reg_stage_dir / vol_id
            specimen_stage_inversion_dir = inv_stage_dir / vol_id

            transform_file = common.getfile_startswith(specimen_stage_reg_dir, ELX_TRANSFORM_NAME)
            parameter_file = common.getfile_startswith(reg_stage_dir, ELX_PARAM_PREFIX)

            # Create the folder to put the specimen inversion parameter files in.
            inv_stage_dir.mkdir(exist_ok=True)

            # Add the stage to the inversion order config (in reverse order), if not already.
            if reg_stage_dir.name not in stages_to_invert['inversion_order']:
                stages_to_invert['inversion_order'].insert(0, reg_stage_dir.name)

            if clobber:
                common.mkdir_force(specimen_stage_inversion_dir)  # Overwrite any inversion file that exist for a single specimen

            # Each registration directory contains a metadata file, which contains the relative path to the fixed volume
            reg_metadata = cfg_load(specimen_stage_reg_dir / common.INDV_REG_METADATA)
            fixed_volume = (specimen_stage_reg_dir / reg_metadata['fixed_vol']).resolve()

            # Invert the Transform parameters with options for normal image inversion

            job = {
                'specimen_stage_inversion_dir': specimen_stage_inversion_dir,
                'parameter_file': abspath(parameter_file),
                'transform_file': transform_file,
                'fixed_volume': fixed_volume,
                'param_file_output_name': 'inversion_parameters.txt',
                'image_replacements': IMAGE_REPLACEMENTS,
                'label_replacements': LABEL_REPLACEMENTS,
                'image_transform_file': IMAGE_INVERTED_TRANSFORM,
                'label_transform_file': LABEL_INVERTED_TRANFORM,
                'clobber': clobber,
                'threads': threads
            }

            jobs.append(job)

    # Run the inversion jobs. Currently using only one thread as it seems that elastix now uses multiple threads on the
    # Inversions

    logging.info('inverting with {} threads: '.format(threads))
    pool = Pool(1) # 17/09/18 If we can get multithreded inversion in elastix 4.9 we can remove the python multithreading
    try:
        pool.map(_invert_transform_parameters, jobs)

    except KeyboardInterrupt:
        print('terminating inversion')
        pool.terminate()
        pool.join()

    # TODO: Should we replace the need for this invert.yaml?
    reg_dir = Path(os.path.relpath(reg_stage_dir, inv_outdir))
    stages_to_invert['registration_directory'] = str(reg_dir)  # Doc why we need this
    # Create a yaml config file so that inversions can be run seperatley
    invert_config = config['inverted_transforms'] / INVERT_CONFIG

    with open(invert_config, 'w') as yf:
        yf.write(yaml.dump(dict(stages_to_invert), default_flow_style=False))


def _invert_transform_parameters(args: Dict):
    """
    Generate a single inverted elastix transform parameter file. This can then be used to invert labels, masks etc.
    If any of the step fail, return as subsequent steps will also fail. The logging of failures is handled
    within each function
    """

    # If we have both the image and label inverted transforms, don't do anything if noclobber is True
    clobber = args['clobber']
    threads = args['threads']

    image_transform_param_path = abspath(join(args['specimen_stage_inversion_dir'], args['image_transform_file']))
    label_transform_param_path = abspath(join(args['specimen_stage_inversion_dir'], args['label_transform_file']))

    if not clobber and isfile(label_transform_param_path) and isfile(image_transform_param_path):
        logging.info('skipping {} as noclobber is True and inverted parameter files exist')
        return

    # Modify the elastix registration input parameter file to enable inversion (Change metric and don't write image results)
    inversion_params = abspath(join(args['specimen_stage_inversion_dir'], args['param_file_output_name'])) # The elastix registration parameters used for inversion
    make_elastix_inversion_parameter_file(abspath(args['parameter_file']), inversion_params, args['image_replacements'])  # I don't think we need the replacements here!!!!!!!!

     # Do the inversion, making the inverted TransformParameters file
    fixed_vol = args['fixed_volume']
    forward_tform_file = abspath(args['transform_file'])
    invert_param_dir = args['specimen_stage_inversion_dir']

    if not invert_elastix_transform_parameters(fixed_vol, forward_tform_file, inversion_params, invert_param_dir, threads):
        return

    # Get the resulting TransformParameters file, and create a transform file suitable for inverting normal volumes
    image_inverted_tform = abspath(join(args['specimen_stage_inversion_dir'], 'TransformParameters.0.txt'))

    if not _modify_inverted_tform_file(image_inverted_tform, image_transform_param_path):
        return

    # Get the resulting TransformParameters file, and create a transform file suitable for inverting label volumes

    # replace the parameter in the image file with label-specific parameters and save in new file. No need to
    # generate one from scratch
    if not make_elastix_inversion_parameter_file(image_transform_param_path, label_transform_param_path, args['label_replacements']):
        return

    _modify_inverted_tform_file(label_transform_param_path)


def get_reg_dirs(config: LamaConfig) -> List[Path]:
    """
    Get the registration output directories paths in the order they were made
    """

    reg_stages = []

    for i, reg_stage in enumerate(config['registration_stage_params']):
        stage_id = reg_stage['stage_id']
        stage_dir = config['root_reg_dir'] / stage_id
        reg_stages.append(stage_dir)
    return reg_stages


def make_elastix_inversion_parameter_file(elx_param_file: Path, newfile_name: str, replacements: Dict):
    """
    Modifies the elastix input parameter file that was used in the original transformation.
    Adds DisplacementMagnitudePenalty (which is needed for inverting)
    Turns off writing the image results at the end as we only need an inverted output file.
    Also changes interpolation order in the case of inverting labels

    Parameters
    ----------
    elx_param_file: str
        path to elastix input parameter file
    newfile_name: str
        path to save modified parameter file to

    """

    try:
        with open(elx_param_file) as old, open(newfile_name, "w") as new:

            for line in old:
                if line.startswith("(Metric "):
                    line = '(Metric "DisplacementMagnitudePenalty")\n'
                if line.startswith('(WriteResultImage '):
                    line = '(WriteResultImage "false")\n'
                if line.startswith('WriteResultImageAfterEachResolution '):
                   continue
                try:
                    param_name = line.split()[0][1:]
                except IndexError:
                    continue  # comment?

                if param_name in replacements:
                    value = replacements[param_name]
                    try:
                        int(value)
                    except ValueError:
                        # Not an int, neeed quotes
                        line = '({} "{}")\n'.format(param_name, value)
                    else:
                        # An int, no quotes
                        line = '({} {})\n'.format(param_name, value)
                new.write(line)
    except IOError as e:
        logging.error("Error modifying the elastix parameter file: {}".format(e))
        return False
    return True


def invert_elastix_transform_parameters(fixed: Path, tform_file: Path, param: Path, outdir: Path, threads: str):
    """
    Invert the transform and get a new transform file
    """
    if not common.test_installation('elastix'):
        raise OSError('elastix not installed')

    cmd = ['elastix',
           '-t0', tform_file,
           '-p', param,
           '-f', fixed,
           '-m', fixed,
           '-out', outdir,
           '-threads', threads   # 11/09/18. This was set to 1. Can iversions take advantage of multithreading?
           ]


    try:
        subprocess.check_output(cmd)
    except (Exception, subprocess.CalledProcessError) as e:
        msg = f'Inverting transform file failed. cmd: {cmd}\n{str(e)}:'
        logging.error(msg)
        logging.exception(msg)
        return False
    return True


def _modify_inverted_tform_file(elx_tform_file: Path, newfile_name: str=None):
    """
    Remove "NoInitialTransform" from the output transform parameter file
    Set output image format to unsigned char. Writes out a modified elastix transform parameter file
    that can be used for inverting volumes

    Parameters
    ----------
    elx_tform_file: str
        path to elastix transform file
    newfile_mame: str
        path to save modified transform file
    """

    if not newfile_name:  # Write to temporary file before overwriting
        new_file = tempfile.NamedTemporaryFile().name
    else:
        new_file = newfile_name

    try:

        with open(new_file, "w+") as new_tform_param_fh,  open(elx_tform_file, "r") as tform_param_fh:

            for line in tform_param_fh:
                if line.startswith('(InitialTransformParametersFileName'):
                    line = '(InitialTransformParametersFileName "NoInitialTransform")\n'
                new_tform_param_fh.write(line)
            new_tform_param_fh.close()
            tform_param_fh.close()

    except IOError:
        logging.warning("Error reading or writing transform files {}".format(elx_tform_file))
        return False

    return True


# def is_euler_stage(tform_param):
#     """
#     Return True if the registration used to create this param file was a Euler transform. Can't currently invert
#     Euler transforms with this method, and is usually not required
#     :param tform_param:
#     :return:
#     """
#     with open(tform_param, 'r') as fh:
#         line = fh.readline()
#         if 'EulerTransform' in line:
#             return True
#         else:
#             return False