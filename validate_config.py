from os.path import join
import os
import sys
from . import common
import logging
import numpy as np
import difflib
from sys import version_info

KNOWN_OPTIONS = (
    'no_qc', 'pad_dims', 'threads', 'filetype', 'compress_averages', 'fixed_volume',  'voxel_size', 'output_dir',
    'generate_new_target_each_stage', 'skip_transform_inversion',  'global_elastix_params', 'registration_stage_params',
    'fixed_mask', 'stats_mask', 'pairwise_registration', 'isosurface_dir', 'label_map', 'inverted_isosurfaces',
    'restart_at_stage', 'label_names', 'generate_deformation_fields', 'inputs', 'skip_deformation_fields',
    'normalisation_roi', 'staging', 'staging_volume', 'histogram_normalise_target', 'data_type', 'glcm',
    'littermate_pattern', 'use_auto_staging'
)

DATA_TYPE_OPTIONS = ('uint8', 'int8', 'int16', 'uint16', 'float32')


def validate_reg_config(config, config_dir):
    """
    Do some checks on the config file to check for errors.

    Parameters
    ----------
    config:

    TODO: Add stats checking. eg. if using lmR, need to specify formula
    """

    # Check if there are any unkown options in the config in order to spot typos
    unkown_options = check_for_unkown_options(config)
    if unkown_options:
        param, suggestions = unkown_options
        if not suggestions:
            suggestions = ["?"]
        logging.error("The following option is not recognised: {}\nDid you mean: {} ?".format(param, ", ".join(suggestions)))
        sys.exit()

    convert_image_pyramid(config)
    pairwise_check(config)

    required_params = ['global_elastix_params',
                       'registration_stage_params']

    if not config.get('pairwise_registration'):
        required_params.append('fixed_volume')

    if config.get('staging'):
        check_staging(config)

    required_present = True
    for p in required_params:
        if p not in config:
            required_present = False
            logging.error("Entry '{}' is required in the config file".format(p))
    if not required_present:
        sys.exit(1)

    validate_filetype(config)

    # If label path or isosurfaces inputs are set, add default out dirs if not set
    if config.get('isosurface_dir'):
        if not config.get('inverted_isosurfaces'):
            config['inverted_isosurfaces'] = 'inverted_isosurfaces'

    if config.get('label_map'):
        if not config.get('inverted_labels'):
            config['inverted_labels'] = 'inverted_labels'

    # Check paths
    if not config.get('inputs'):
        config['inputs'] = join(config_dir, 'inputs')
    else:
        config['inputs'] = join(config_dir, config['inputs'])

    paths = []
    if not config.get('restart_at_stage'):  # if not,  we do not need inputs
        paths.append(config.get('inputs'))

    if not config.get('pairwise_registration'):
        paths.append(config.get('fixed_volume'))

    if config.get('fixed_mask'):
        paths.append(config.get('fixed_mask'))
    if config.get('stats_mask'):
        paths.append(config.get('stats_mask'))
    if config.get('label_map'):
        paths.append(config.get('label_map'))
    if config.get('organ_names'):
        paths.append(config.get('organ_names'))
    if config.get('isosurface_dir'):
        paths.append(config.get('isosurface_dir'))

    failed_paths = check_paths(config_dir, paths)
    if len(failed_paths) > 0:
        for f in failed_paths:
            logging.error("Cannot find '{}'. All paths need to be relative to config file".format(f))
        sys.exit(1)

    stages = config['registration_stage_params']

    # Check for correct deformation fields generation parameter
    if config.get('generate_deformation_fields'):
        for _, def_info in config.get('generate_deformation_fields').items():
            stage_id = def_info[0]
            def_stage_found = False
            for stage in stages:
                if stage['stage_id'] == stage_id:
                    def_stage_found = True
                    # Check the resolutions make sense
                    # if len(def_info) > 1:
                    #     if not config['global_elastix_params'].get('WriteTransformParametersEachResolution'):
                    #         logging.error("If specifiying resolutions to generate deformations from, 'WriteTransformParametersEachResolution: true' must be added to 'global_elastix_params' ")
                    #         sys.exit(1)
                    #     for res in def_info[1:]:
                    #         try:
                    #             int(res)
                    #         except ValueError:
                    #             logging.error("Resolutions in 'generate_deformation_fields' should be integers")
                    #             sys.exit(1)
                    #     final_res_def = max( def_info[1:])
                    #     final_res_reg = stage['elastix_parameters']['NumberOfResolutions']
                    #     if final_res_def > final_res_reg:
                    #         logging.error("There are only {} resolutions specified in {}. But the final resolution specified in 'generate_deformation_fields' is {} ".format(final_res_reg, stage_id, final_res_def))
                    #         sys.exit(1)

            if not def_stage_found:
                logging.error("'generate_deformation_fields' should refer to a registration 'stage_id'")
                sys.exit(1)

    # check we have some registration stages specified
    if len(stages) < 1:
        logging.error("No registration stages specified")
        sys.exit(1)

    # Check normalised output ROI
    normalised_roi = config.get('normalisation_roi')
    if normalised_roi:
        if normalised_roi != 'mask':
            starts = normalised_roi[0]
            ends = normalised_roi[1]
            if starts[0] >= ends[0] or starts[1] >= ends[1] or starts[2] >= ends[2]:
                logging.error("Error in normalization ROI. should be [x start, y start , z start], [x end, y end, z end] or 'mask'")
                sys.exit(1)
    # Some stuff in the satges
    num_def_stages = 0
    for stage in stages:
        inherit_id = stage.get('inherit_elx_params')
        if inherit_id:
            found_id = False
            for s in stages:
                if s.get('stage_id') == inherit_id:
                    found_id = True
                    break
            if not found_id:
                logging.error("Could not find the registration stage to inherit from '{}'".format(inherit_id))
                sys.exit(1)

    voxel_size = config.get('voxel_size')
    if voxel_size:
        try:
            float(voxel_size)
        except ValueError:
            logging.error("'voxel_size' should be a number")
            sys.exit(1)

    pad_dims = config.get('pad_dims')
    if not isinstance(pad_dims, bool):
        if not isinstance(pad_dims, str):
            if len(pad_dims) != 3:
                logging.error('Pad dims should be either true, false, or a list of x,y,z dimensions. eg [100, 200, 300]')
                sys.exit(1)
        else:
            logging.error('Pad dims should be either true, false, or a list of x,y,z dimensions. eg [100, 200, 300]')
            sys.exit(1)

    check_images(config_dir, config)


def check_images(config_dir, config):
    # Check for image paths
    img_dir = join(config_dir, config['inputs'])
    # Inputs is a folder
    if os.path.isdir(img_dir):
        imgs = os.listdir(img_dir)
    # Inputs is a list of paths
    elif os.path.isfile(img_dir):
        imgs = common.get_inputs_from_file_list(img_dir, config_dir)
    else:
        logging.error("'inputs:' should refer to a sirectory of images or a file containing image paths")
        sys.exit(1)
    logging.info('validating input volumes')

    if not config.get('restart_at_stage'):  # No need to validate volumes, they were created by lama
        dtypes = {}

        for im_name in imgs:
            image_path = join(img_dir, im_name)

            array_load = common.LoadImage(image_path)
            if not array_load:
                logging.error(array_load.error_msg)
                sys.exit(1)

            check_dtype(config, array_load.array, array_load.img_path)
            check_16bit_elastix_parameters_set(config, array_load.array)
            dtypes[im_name] = array_load.array.dtype

        if len(set(dtypes.values())) > 1:
            dtype_str = ""
            for k, v in list(dtypes.items()):
                dtype_str += k + ':\t' + str(v) + '\n'
            logging.warning('The input images have a mixture of data types\n{}'.format(dtype_str))


def check_16bit_elastix_parameters_set(config, array):
    if array.dtype in (np.int16, np.uint16):
        try:
            internal_fixed = config['global_elastix_params']['FixedInternalImagePixelType']
            internal_mov = config['global_elastix_params']['MovingInternalImagePixelType']
            if internal_fixed != 'float' or internal_mov != 'float':
                raise TypeError
        except (TypeError, KeyError):
            logging.warning(
                "If using 16 bit input volumes, 'FixedInternalImagePixelType' and 'MovingInternalImagePixelType should'" \
                "be set to 'float' in the global_elastix_params secion of the config file")
            sys.exit(1)


def check_dtype(config, array, img_path):

    # Check that bit depth is correct
    # TODO fix the bit deph validation
    data_type = config.get('data_type')
    if data_type:  # Currently only checking for int8 and int16. Add float checking as well

        if data_type not in DATA_TYPE_OPTIONS:
            sys.exit('Data type must be one of {}'
                     '\nleave out data_type or set to None if data_type checking not required'.
                     format(str(DATA_TYPE_OPTIONS)))

        if not np.issubdtype(data_type, array.dtype):
            raise ValueError('data type given in config is:{}\nThe datatype for image {} is {}'.
                             format(data_type, img_path, array.dtype))


def check_paths(config_dir, paths):
    failed = []
    for p in paths:
        path = os.path.join(config_dir, p)
        if not os.path.exists(path):
            failed.append(p)
    return failed


def check_staging(config):
    """
    TODO need to add check for mask and label for 2 methods
    Parameters
    ----------
    config

    Returns
    -------

    """
    stage_types = list(common.CONFIG_OPPS['staging'].values())

    if config.get('staging') not in stage_types:
        sys.exit("{} is not a valid staging method".format(config['staging'].get('method')))

    if config.get('staging') == 'label_length' and config.get("skip_transform_inversion"):
        sys.exit("For label length staging to work, labels need to be iverted."
                 " Ensure 'skip_transform_inversion' is not set to 'true'")


def check_for_unkown_options(config):
    for param in config:
        if param not in KNOWN_OPTIONS:
            closest_match = difflib.get_close_matches(param, KNOWN_OPTIONS)
            return param, closest_match
    return False


def validate_filetype(config):
    """
    Filetype can be sepcified in the elastix config section, but this intereferes with LAMA config section
    Parameters
    ----------
    config: dict
        The main config
    -------
    """
    if not config.get('filetype'):
        config['filetype'] = 'nrrd'
    if config['filetype'] not in ('nrrd', 'nii'):
        logging.error("'filetype' should be 'nrrd' or 'nii")
        sys.exit(1)
    gep = config.get('global_elastix_params')
    gep['ResultImageFormat'] = config['filetype']


def pairwise_check(config):
    """
    Do not generate images for pairwise registration as the images are generated by applying mean transform
    Parameters
    ----------
    config
    -------

    """
    if config.get('pairwise_registration'):
        gep = config.get('global_elastix_params')
        gep['WriteResultImage'] = 'false'
        gep['WriteResultImageAfterEachResolution'] = 'false'


def convert_image_pyramid(config):
    """
    The elastix image pyramid needs to be specified for each dimension for each resolution

    This function allows it to specified for just one resolution as it should always be the same for each dimension.
    Convert to elastix required format
    Parameters
    ----------
    config: dict
        paramters

    """
    for stage in config['registration_stage_params']:
        elx_params = stage['elastix_parameters']
        if elx_params.get('ImagePyramidSchedule') and elx_params.get('NumberOfResolutions'):
            num_res = int(elx_params.get('NumberOfResolutions'))
            lama_schedule = elx_params.get('ImagePyramidSchedule')
            if len(lama_schedule) == num_res:
                elastix_shedule = []
                for i in lama_schedule:
                    elastix_shedule.extend([i, i, i])
                elx_params['ImagePyramidSchedule'] = elastix_shedule

