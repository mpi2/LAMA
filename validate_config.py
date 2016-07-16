from os.path import join
import os
import sys
import common
import logging
import numpy as np


def validate_reg_config(config, config_dir):
    """
    Do some checks on the config file to check for errors.
    :param config:
    :return: list, empty or containing errors


    TODO: Add stats checking. eg. if using lmR, need to specify formula
    """
    report = []
    required_params = ['global_elastix_params',
                       'registration_stage_params']

    if not config.get('pairwise_registration'):
        required_params.append('fixed_volume')

    for p in required_params:
        if p not in config:
            report.append("Entry '{}' is required in the config file".format(p))

    if not config.get('filetype'):
        config['filetype'] = 'nrrd'

    # If label path or isosurfaces inputs are set, add default out dirs if not set
    if config.get('isosurface_dir'):
        if not config.get('inverted_isosurfaces'):
            config['inverted_isosurfaces'] = 'inverted_isosurfaces'

    if config.get('label_map_path'):
        if not config.get('inverted_labels'):
            config['inverted_labels'] = 'inverted_labels'

    # Check paths
    if not config.get('inputvolumes_dir'):
        config['inputvolumes_dir'] = join(config_dir, 'inputs')
    paths = [config.get('inputvolumes_dir')]

    if not config.get('pairwise_registration'):
        paths.append(config.get('fixed_volume'))

    if config.get('fixed_mask'):
        paths.append(config.get('fixed_mask'))
    if config.get('label_map_path'):
        paths.append(config.get('label_map_path'))
    if config.get('organ_names'):
        paths.append(config.get('organ_names'))
    if config.get('isosurface_dir'):
        paths.append(config.get('isosurface_dir'))

    failed_paths = check_paths(config_dir, paths)
    if len(failed_paths) > 0:
        for f in failed_paths:
            report.append("Cannot find '{}'. All paths need to be relative to config file".format(f))

    stages = config['registration_stage_params']

    # Check for correct deformation fields generation parameter
    # TODO: rewrite for the new seformation config specificaiton
    # def_start_stage = config.get('generate_deformation_fields')
    # if def_start_stage:
    #     def_stage_found = False
    #     for stage in stages:
    #         if stage['stage_id'] == def_start_stage:
    #             def_stage_found = True
    #     if not def_stage_found:
    #         report.append("Error: 'generate_deformation_fields' should refer to a registration 'stage_id'")

    # check we have some registration stages specified
    if len(stages) < 1:
        report.append("No stages specified")

    # Check normalised output ROI
    normalised_roi = config.get('background_roi_zyx_norm')
    if normalised_roi:
        starts = normalised_roi[0]
        ends = normalised_roi[1]
        if starts[0] >= ends[0] or starts[1] >= ends[1] or starts[2] >= ends[2]:
            report.append("Error in normalization ROI. should be [z start, y start , x start], [z end, y end, x end]")

    # Some stuff in the satges
    num_def_stages = 0
    for stage in stages:
        if stage.get('generate_deformation_fields'):
            num_def_stages += 1
        inherit_id = stage.get('inherit_elx_params')
        if inherit_id:
            found_id = False
            for s in stages:
                if s.get('stage_id') == inherit_id:
                    found_id = True
                    break
            if not found_id:
                report.append("Could not find the stage to inherit from '{}'".format(inherit_id))

    if num_def_stages > 1:
        report.append("multiple 'generate_deformation_fields' entries found. There should only be one at a certain registration stage ")

    # Check whether images are 16 bit and if so whether internal representation is set to float
    img_dir = join(config_dir, config['inputvolumes_dir'])
    if os.path.isdir(img_dir):
        imgs = os.listdir(img_dir)
    else:
        imgs = common.get_inputs_from_file_list(img_dir, config_dir)

    logging.info('validating input volumes')
    if not config.get('restart_at_stage'):
        for im_name in imgs:
            image_path = join(img_dir, im_name)

            array_load = common.LoadImage(image_path)
            if not array_load:
                logging.error(array_load.error_msg)
                sys.exit()
            if array_load.array.dtype in (np.int16, np.uint16):
                try:
                    internal_fixed = config['global_elastix_params']['FixedInternalImagePixelType']
                    internal_mov = config['global_elastix_params']['MovingInternalImagePixelType']
                    if internal_fixed != 'float' or internal_mov != 'float':
                        raise TypeError
                except (TypeError, KeyError):
                    report.append("If using 16 bit input volumes, 'FixedInternalImagePixelType' and 'MovingInternalImagePixelType should'" \
                                  "be set to 'float' in the global_elastix_params secion of the config file")

    check_non_path_options(config, report)

    if len(report) > 0:
        for r in report:
            logging.error(r)
        sys.exit('Registration config file error. PLease see LAMA.log')
    return True

def check_non_path_options(config, report):
    voxel_size = config.get('voxel_size')
    if voxel_size:
        try:
            float(voxel_size)
        except ValueError:
            report.append('voxel_size should be a number')

    pad_dims = config.get('pad_dims')
    if not isinstance(pad_dims, bool):
        if not isinstance(pad_dims, basestring):
            if len(pad_dims) != 3:
                report.append('Pad dims should be either true, false, or a list of x,y,z dimensions. eg [100, 200, 300]')
        else:
            report.append('Pad dims should be either true, false, or a list of x,y,z dimensions. eg [100, 200, 300]')



def check_paths(config_dir, paths):
    failed = []
    for p in paths:
        path = os.path.join(config_dir, p)
        if not os.path.exists(path):
            failed.append(p)
    return failed
