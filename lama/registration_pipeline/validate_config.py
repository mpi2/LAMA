from os.path import join
import os
import sys
from lama import common
from lama import paths
from logzero import logger as logging
import numpy as np
import difflib
from pathlib import Path
from collections import OrderedDict
from typing import Union, List, Dict
from sys import version_info
import toml

from lama.staging.staging_metric_maker import STAGING_METHODS, DEFAULT_STAGING_METHOD

outdir_path_names = ['isosurface_dir', 'inverted_isosurfaces']

DATA_TYPE_OPTIONS = ('uint8', 'int8', 'int16', 'uint16', 'float32')


class LamaConfig:
    """
    Contains information derived form the user config such as paths
    """

    def __init__(self, config_path):

        # read in the the config
        with open(config_path) as fh:
            self.config = toml.load(fh)

        self.required_keys = ['global_elastix_params',
                         'registration_stage_params']

        self.output_path_names = OrderedDict({
            'target_folder': 'target',
            'output_dir': 'output',
            'deformations': 'deformations',
            'jacobians': 'jacobians',
            'jacmat': 'jacobian_matrices',
            'glcms': 'glcms',
            'root_reg_dir': 'registrations',
            'inverted_labels': 'inverted_labels',
            'inverted_stats_masks': 'inverted_stats_masks',
            'organ_vol_result_csv': common.ORGAN_VOLUME_CSV_FILE
        })

        self.target_names = (
            'fixed_mask',
            'stats_mask',
            'fixed_volume',
            'label_map',
            'label_names'
        )

        self.input_options = {

            # parameter: [options...], default]
            # Opions can be types or functions
            'no_qc': (bool, False),
            'threads': (int, 4),
            'filetype': self.validate_filetype,
            'voxel_size': (float, 14.0),
            'generate_new_target_each_stage': (bool, False),
            'skip_transform_inversion': (bool, False),
            'pairwise_registration': (bool, False),
            'generate_deformation_fields': (dict, None),
            'skip_deformation_fields': (bool, True),
            'staging': (list(STAGING_METHODS.keys()), DEFAULT_STAGING_METHOD),
            'data_type': (['uint8', 'int8', 'int16', 'uint16', 'float32'], 'uint8'),
            'glcm': (bool, False),
            'config_version': (float, None)
        }

        self.all_keys = list(self.required_keys) + list(self.output_path_names.keys()) + list(self.target_names) + list(self.input_options.keys())

        # options is where the final options (either default or from config) are stored.
        # Paths from config or default will have been resolved relative to config directoryu
        self.options = {}

        self.config_dir = config_path.parent

        # Check if there are any unkown options in the config in order to spot typos
        self.check_for_unkown_options()

        self.check_required()

        self.convert_image_pyramid()

        self.pairwise_check()

        self.check_paths()

        self.check_options()

        self.check_stages()

        self.check_images()

        self.resolve_output_paths()

    def __getitem__(self, item):
        return self.option

    def resolve_output_paths(self):
        for k, p in self.output_path_names.items():

            if k in ('output_dir', 'target_folder'):
                path = self.config_dir  / p
            else:
                path = self.config_dir / self.options['output_dir'] / p

            self.options[k] = path

    def check_options(self):
        """
        Check the options in the lama section of the config. Perform apppropriate validation, or set default
        add the options to self.options

        """
        for option, validation in self.input_options.items():

            # Run a function to check
            if callable(validation):
                validation()  # Should have it's own merror message and exit

            else:
                value = self.config.get(option, None)

                checker = validation[0]
                default = validation[1]

                if value is None:
                    value = default
                    continue

            # Check bool
            if checker == bool:
                if type(value) != bool:
                    sys.exit(f'{option} should be a bool not a {type(value)}')

            # validate by type
            elif type(checker) == type:
                if not isinstance(value, checker):
                    sys.exit(f'"{option}" should be type {str(checker)}')


            # Check for a list of options
            else:
                try:
                    iter(checker)
                except TypeError:
                    pass
                else:
                    if value not in checker:
                        sys.exit(f'{option} should be one of {checker}')

        self.options[option] = value

    def path_mkdir(self) -> Path:
        """
        Get a a path and make the directory as well including parents

        Returns
        -------

        """
        pass

    # Volumes that are in the population average (target) space whose names should be specified using the following keys

    # Paths in the output directory have the following defaults, but can be overridden in the config file


    def validate_filetype(self):
        """
        Filetype can be sepcified in the elastix config section, but this intereferes with LAMA config section
        Parameters
        ----------
        config: dict
            The main config
        -------
        """
        # todo: this is a bit broke, are gloval_elx.. in optionns or config at this point
        ft = self.config.get('filetype', 'nrrd')

        if ft not in ('nrrd', 'nii'):
            logging.error("'filetype' should be 'nrrd' or 'nii")
            sys.exit(1)

        self.options['filetype'] = ft

        gep = self.config.get('global_elastix_params')
        gep['ResultImageFormat'] = ft

    def check_stages(self):
        config = self.config

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


    def check_required(self):
        if not self.config.get('pairwise_registration'):
            # No fixed volume if we are doing pairwise registration
            self.required_keys.append('fixed_volume')

        failed = False
        for k in self.required_keys:
            if k not in self.config:
                failed = True
                logging.error("Entry '{}' is required in the config file".format(k))
        if failed:
            sys.exit(1)
        else:
            self.options[k] = self.config[k]

    def check_images(self):
        """
        validate that image paths are correct and give loadeable volumes
        """

        img_dir = self.options['inputs']

        # Inputs is a folder
        if os.path.isdir(img_dir):
            imgs = os.listdir(img_dir)

        # Inputs is a list of paths
        elif os.path.isfile(img_dir):
            imgs = common.get_inputs_from_file_list(img_dir, self.config_dir)
        else:
            logging.error("'inputs:' should refer to a sirectory of images or a file containing image paths")
            sys.exit(1)
        logging.info('validating input volumes')

        dtypes = {}

        for im_name in imgs:
            image_path = join(img_dir, im_name)

            array_load = common.LoadImage(image_path)
            if not array_load:
                logging.error(array_load.error_msg)
                raise FileNotFoundError(f'cannot load {image_path}')

            self.check_dtype(self.config, array_load.array, array_load.img_path)
            self.check_16bit_elastix_parameters_set(self.config, array_load.array)
            dtypes[im_name] = array_load.array.dtype

        if len(set(dtypes.values())) > 1:
            dtype_str = ""
            for k, v in list(dtypes.items()):
                dtype_str += k + ':\t' + str(v) + '\n'
            logging.warning('The input images have a mixture of data types\n{}'.format(dtype_str))


    def check_16bit_elastix_parameters_set(self, config, array):
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

    def check_dtype(self, config, array, img_path):

        # Check that bit depth is correct
        # TODO fix the bit deph validation
        data_type = config.get('data_type')
        if data_type:  # Currently only checking for int8 and int16. Add float checking as well

            if data_type not in DATA_TYPE_OPTIONS:
                sys.exit('Data type must be one of {}'
                         '\nleave out data_type or set to None if data_type checking not required'.
                         format(str(DATA_TYPE_OPTIONS)))

            if not np.issubdtype(data_type, array.dtype):
                msg = "data type given in config is:{}\nThe datatype for image {} is {}".format(data_type, img_path, array.dtype)
                logging.error(msg)
                raise ValueError(msg)

    def check_paths(self):
        """
        Check the presence of files named in the config.
        All paths are relative to the config dir. Set the resolved paths in self.options
        """
        failed = []

        # Check the target paths first. Paths that should be in the target folder
        target_folder = self.config_dir / self.config['target_folder']
        if not target_folder.is_dir():
            sys.exit(f'target foldet not found: {target_folder}')

        self.options['target_dir'] = target_folder

        for tn in self.target_names:
            if tn in self.config:
                target_path = target_folder / self.config[tn]
                if not target_path.exists():
                    failed.append(target_path)
                else:
                    self.options[tn] = target_path
            else:
                self.options[tn] = None

        # Check paths
        if not self.config.get('inputs'):
            #
            self.options['inputs'] = self.config_dir / 'inputs'
        else:
            self.options['inputs'] = self.config_dir / self.config['inputs']

        if failed:
            for f in failed:
                logging.error("Cannot find '{}'. All paths need to be relative to config file".format(f))
            raise ValueError("Cannot find '{}'. All paths need to be relative to config file".format(f))

    def check_for_unkown_options(self):

        for param in self.config:
            if param not in self.all_keys:
                closest_match = difflib.get_close_matches(param, self.all_keys)

                if not closest_match:
                    suggestions = ["?"]
                logging.error(
                    "The following option is not recognised: {}\nDid you mean: {} ?".format(param, ", ".join(suggestions)))
                sys.exit()


    def pairwise_check(self):
        """
        Do not generate images for pairwise registration as the images are generated by applying mean transform
        Parameters
        ----------
        config
        -------

        """
        if self.config.get('pairwise_registration'):
            gep = self.config.get('global_elastix_params')
            gep['WriteResultImage'] = 'false'
            gep['WriteResultImageAfterEachResolution'] = 'false'


    def convert_image_pyramid(self):
        """
        The elastix image pyramid needs to be specified for each dimension for each resolution

        This function allows it to specified for just one resolution as it should always be the same for each dimension.
        Convert to elastix required format

        Parameters
        ----------
        config: dict
            paramters

        """
        for stage in self.config['registration_stage_params']:
            elx_params = stage['elastix_parameters']
            if elx_params.get('ImagePyramidSchedule') and elx_params.get('NumberOfResolutions'):
                num_res = int(elx_params.get('NumberOfResolutions'))
                lama_schedule = elx_params.get('ImagePyramidSchedule')
                if len(lama_schedule) == num_res:
                    elastix_shedule = []
                    for i in lama_schedule:
                        elastix_shedule.extend([i, i, i])
                    elx_params['ImagePyramidSchedule'] = elastix_shedule

