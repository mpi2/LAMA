from os.path import join
import os
import sys
import shutil
import difflib
from pathlib import Path
from collections import OrderedDict

from logzero import logger as logging
import numpy as np
import toml

from lama import common
from lama.staging.staging_metric_maker import STAGING_METHODS, DEFAULT_STAGING_METHOD


DATA_TYPE_OPTIONS = ('uint8', 'int8', 'int16', 'uint16', 'float32')


class LamaConfigError(BaseException):

    pass


class LamaConfig:
    """
    Contains information derived form the user config such as paths
    The keys from output_paths, target_names and input_options will be avaialble via the __getitem__
    """

    def __init__(self, config_path):
        """
        Parameters
        ----------
        config_path
            pat to the lama config file

        Raises
        ------
        OSError of subclasses thereof if config file cannot be opened
        """

        self.config_path = config_path
        # read in the the config
        with open(config_path) as fh:
            self.config = toml.load(fh)


        # The variable names mapped to the actual names of output directories
        # If the value is a string, it will be created in the output_dir
        # If the value is a tuple [0] is the folder name and the rest are parent folders
        self.output_path_names = OrderedDict({
            # output_dir must always be 'output' as some other modules depend upon this
            # Must add way to enforce as it can be overriden in the config at the moment
            'output_dir': 'output',
            'target_folder': 'target',
            'qc_dir': 'qc',
            'input_image_histograms': ('input_image_histograms', 'qc'),
            'qc_registered_images': ('qc_registered_images', 'qc'),
            'metric_charts_dir': ('metric_charts', 'qc'),
            'registered_midslice_dir': ('registered_midslices', 'qc'),
            'inverted_label_overlay_dir': ('inverted_label_overlay', 'qc'),
            'cyan_red_dir': ('cyan_red_overlay', 'qc'),
            'average_folder': 'averages',
            'deformations': 'deformations',
            'jacobians': 'jacobians',
            'log_jacobians': 'log_jacobians',
            'jacmat': 'jacobian_matrices',
            'glcm_dir': 'glcms',
            'root_reg_dir': 'registrations',
            'inverted_transforms': 'inverted_transforms',
            'inverted_labels': 'inverted_labels',
            'inverted_stats_masks': 'inverted_stats_masks',
            'organ_vol_result_csv': common.ORGAN_VOLUME_CSV_FILE
        })

        # Options in the config that map to files that can be present in the target folder
        self.target_names = (
            'fixed_mask',
            'stats_mask',
            'fixed_volume',
            'label_map',
            'label_names'
        )

        self.input_options = {

            # parameter: [options...], default]
            # Options can be types or functions
            'global_elastix_params': ('dict', 'required'),
            'registration_stage_params': ('dict', 'required'),
            'no_qc': ('bool', False),
            'threads': ('int', 4),
            'filetype': ('func', self.validate_filetype),
            'voxel_size': ('float', 14.0),
            'generate_new_target_each_stage': ('bool', False),
            'skip_transform_inversion': ('bool', False),
            'pairwise_registration': ('bool', False),
            'generate_deformation_fields': ('dict', None),
            'skip_deformation_fields': ('bool', True),
            'staging': ('func', self.validate_staging),
            'data_type': (['uint8', 'int8', 'int16', 'uint16', 'float32'], 'uint8'),
            'glcm': ('bool', False),
            'config_version': ('float', 1.1)
        }

        # The paths to each stage utput dir: stage_id: Path
        self.stage_dirs = OrderedDict()

        self.all_keys = list(self.output_path_names.keys()) + list(self.target_names) + list(self.input_options.keys())

        # options is where the final options (either default or from config) are stored.
        # Paths from config or default will have been resolved relative to config directoryu
        self.options = {}

        self.config_dir = config_path.parent

        # Check if there are any unkown options in the config in order to spot typos
        self.check_for_unknown_options()

        self.convert_image_pyramid()

        self.pairwise_check()

        self.check_paths()

        self.check_options()

        # self.check_images()

        self.resolve_output_paths()

        self.check_stages()

    def __getitem__(self, item):
        return self.options[item]

    def resolve_output_paths(self):
        """
        Make absolute paths from the self.output_path_names dict.
        Paths will be below config_dir / output_dir

        Add the absolute paths to self.options

        """
        for folder_var, folder_name in self.output_path_names.items():

            if folder_var in ('output_dir', 'target_folder'):
                path = self.config_dir / folder_name
            else:
                if isinstance(folder_name, tuple):
                    name, *parents = folder_name
                    path = (self.config_dir / self.options['output_dir']).joinpath(*parents) / name
                else:
                    path = self.config_dir / self.options['output_dir'] / folder_name

            self.options[folder_var] = path

    def check_options(self):
        """
        Check the options in the lama section of the config. Perform apppropriate validation, or set default
        add the options to self.options

        """
        for option, validation in self.input_options.items():

            # # Run a function to check
            # if callable(validation):
            #     validation()  # Should have it's own merror message and exit

            checker = validation[0]
            default = validation[1]

            if default == 'required':
                if option not in self.config:
                    raise LamaConfigError(f'{option} is a required option')
            value = self.config.get(option, default)

            if checker == 'bool':
                if type(value) != bool:
                    raise (f'{option} should be a bool not a {type(value)}')

            elif checker == 'float':
                if not isinstance(value, float):
                    raise LamaConfigError(f'"{option}" should be a number (float. eg 4.6)')

            elif checker == 'int':
                isint = False
                if isinstance(value, int):
                    isint = True
                elif isinstance(value, float):
                    if value / 1 == 0:  # If it's an int in float form, that OK
                        isint = True
                if not isint:
                    raise LamaConfigError(f'"{option}" should be type {str(checker)}')

            elif checker == 'func':
                validation[1]()
                continue # Option set in function

            # Check for a list of options
            elif isinstance(checker, list):
                if value not in checker:
                    raise LamaConfigError(f'{option} should be one of {checker}')

            self.options[option] = value

    def mkdir(self, name, clobber=True) -> Path:
        """
        Get a a path and make the directory as well including parents

        Returns
        -------

        """

        if name in self.options:
            dir_: Path = self.options[name]
            if clobber:
                if dir_.is_dir():
                    shutil.rmtree(dir_)
            dir_.mkdir(exist_ok=True)
        else:
            raise LamaConfigError(f'{name} is not a specified folder')
        return dir_

    def validate_staging(self):
        st = self.config.get('staging')
        if not st:
            st = DEFAULT_STAGING_METHOD
        if st == 'none':
            st = None
        else:
            if st not in list(STAGING_METHODS.keys()):
                raise LamaConfigError('staging must be one of {}'.format(','.join(list(STAGING_METHODS.keys()))))
            if st == 'embryo_volume' and not self.config.get('stats_mask'):
                raise LamaConfigError("To calculate embryo volume a 'stats_mask' should be added to the config ")

        self.options['staging'] = st


    def validate_filetype(self):
        """
        Filetype can be specified in the elastix config section, but this intereferes with LAMA config section
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
            raise LamaConfigError()

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
                    raise LamaConfigError()

        # check we have some registration stages specified
        if len(stages) < 1:
            logging.error("No registration stages specified")
            raise LamaConfigError()

        for stage in stages:

            # Make a path for the output of this stage
            path = self.options['root_reg_dir'] / stage['stage_id']
            self.stage_dirs[stage['stage_id']] = path

            inherit_id = stage.get('inherit_elx_params')
            if inherit_id:
                found_id = False
                for s in stages:
                    if s.get('stage_id') == inherit_id:
                        found_id = True
                        break
                if not found_id:
                    logging.error("Could not find the registration stage to inherit from '{}'".format(inherit_id))
                    raise LamaConfigError()

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
            logging.error("'inputs:' should refer to a directory of images or a file containing image paths")
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
                raise LamaConfigError(msg)

    def check_paths(self):
        """
        Check the presence of files named in the config.
        All paths are relative to the config dir. Set the resolved paths in self.options
        """
        failed = []

        # Check the target paths first. Paths that should be in the target folder
        target_folder = (self.config_dir / self.config['target_folder']).resolve()
        if not target_folder.is_dir():
            raise LamaConfigError(f'target folder not found: {target_folder}')

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
            raise LamaConfigError("Cannot find '{}'. All paths need to be relative to config file".format(f))

    def check_for_unknown_options(self):

        for param in self.config:
            if param not in self.all_keys:
                closest_matches = difflib.get_close_matches(param, self.all_keys)

                if not closest_matches:
                    closest_matches = ["?"]

                msg = "The following option is not recognised: {}\nDid you mean: {} ?".format(param, ", ".join(closest_matches))
                logging.error(msg)
                raise LamaConfigError(msg)

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

