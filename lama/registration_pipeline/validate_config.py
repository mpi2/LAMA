from os.path import join
import os
import sys
import shutil
import difflib
from pathlib import Path
from collections import OrderedDict
from typing import Union, Dict

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
    Contains information derived form the user config such as paths and options.
    The keys from output_paths, target_names and input_options will be avaialble via the __getitem__ function

    Example
    -------
    lc = LamaConfig(Path('cfg_path'))

    """

    def __init__(self, config: Union[Path, Dict], cfg_path: Path=None, no_validate=False):
        """
        Parameters
        ----------
        config
            path to the lama config file
            or
            config dictionary
        cfg_path
            Used for testing. If we want to pass in a dict rather than a path we with also need a path of the project
            directory, which is normally the cfg parent directory

        Raises
        ------
        OSError or subclasses thereof if config file cannot be opened
        """
        if isinstance(config, dict):
            if cfg_path is None:
                raise ValueError("Please supply a project root path")
            self.config = config
            config_path = cfg_path
        elif isinstance(config, Path):
            self.config = common.cfg_load(config)
            config_path = config
        else:
            raise ValueError("config must me a Path or Dict")
        self.config_path = Path(config_path)

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
            'organ_vol_result_csv': common.ORGAN_VOLUME_CSV_FILE,
            'additional_seg_dir': 'additional_seg'

        })

        # Options in the config that map to files that can be present in the target folder
        self.target_names = (
            'fixed_mask',
            'stats_mask',
            'fixed_volume',
            'label_map',
            'label_info'
        )

        self.input_options = {
            # Config parameters to be validated (non-elastix related parameters)
            # parameter: ([options...], default)
            # Options can be types to check against or functions that retrn True is value is valid
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
            'staging': ('func', self.validate_staging),
            'data_type': (['uint8', 'int8', 'int16', 'uint16', 'float32'], 'uint8'),
            'glcm': ('bool', False),
            'config_version': ('float', 1.1),
            'stage_targets': (Path, False),
            'fix_folding': (bool, False),
            # 'inverse_transform_method': (['invert_transform', 'reverse_registration'], 'invert_transform')
            'label_propagation': (['invert_transform', 'reverse_registration'], 'reverse_registration'),
            'skip_forward_registration': (bool, False),
            'seg_plugin_dir': (Path, None),

            # The following options are used for saving dsk space
            'write_deformation_vectors': (bool, False),
            'delete_inverted_transforms': (bool, False),
            'write_raw_jacobians': (bool, False),
            'write_log_jacobians': (bool, True),
        }

        # The paths to each stage output dir: stage_id: Path
        self.stage_dirs = OrderedDict()

        self.all_keys = list(self.output_path_names.keys()) + list(self.target_names) + list(self.input_options.keys())

        # options is where the final options (either default or from config) are stored.
        # Paths from config or default will have been resolved relative to config directoryu
        self.options = {}

        self.config_dir = config_path.parent

        # Check if there are any unkown options in the config in order to spot typos
        if no_validate:
            return

        self.check_for_unknown_options()

        self.convert_image_pyramid()

        self.pairwise_check()

        self.check_paths()

        self.check_options()

        self.check_images()

        self.resolve_output_paths()

        self.check_stages()


    def __getitem__(self, item):
        return self.options[item]

    def __setitem__(self, key, value):
        # For debugging
        self.options[key] = value

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

            if st == 'embryo_volume':
                if not self.config.get('stats_mask') or self.config.get('skip_transform_inversion'):
                    raise LamaConfigError("To calculate embryo volume the following options must be set\n"
                                      "'stats_mask' which is tight mask use for statistical analysis and calcualting whoel embryo volume\n"
                                      "'skip_transform_inversion' must not be False the inversions are needed to calculate embryo volume")

                non_def_stages = [x for x in self.config['registration_stage_params'] if
                                                        x['elastix_parameters'].get('Transform') and
                                                        x['elastix_parameters']['Transform'] in
                                                        ['SimilarityTransform', 'AffineTransform']]
                if not non_def_stages:
                    raise LamaConfigError("In order to calculate embryo volume an affine or similarity stage is needed")

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

        For now just check for spaces in filenames
        """

        img_dir = self.options['inputs']

        # Inputs is a folder
        if os.path.isdir(img_dir):
            imgs = os.listdir(img_dir)

        # Inputs is a list of paths
        elif os.path.isfile(img_dir):
            imgs = common.get_inputs_from_file_list(img_dir, self.config_dir)
        else:
            logging.error("'inputs' directory should be in same directory as config and should contain images or a file "
                          "containing image paths")
            sys.exit(1)
        logging.info('validating input volumes')

        dtypes = {}

        for im_name in imgs:
            image_path = join(img_dir, im_name)

            self.space_filename_check(image_path)

        #     array_load = common.LoadImage(image_path)
        #     if not array_load:
        #         logging.error(array_load.error_msg)
        #         raise FileNotFoundError(f'cannot load {image_path}')
        #
        #     self.check_dtype(self.config, array_load.array, array_load.img_path)
        #     self.check_16bit_elastix_parameters_set(self.config, array_load.array)
        #     dtypes[im_name] = array_load.array.dtype
        #
        # if len(set(dtypes.values())) > 1:
        #     dtype_str = ""
        #     for k, v in list(dtypes.items()):
        #         dtype_str += k + ':\t' + str(v) + '\n'
        #     logging.warning('The input images have a mixture of data types\n{}'.format(dtype_str))

    def space_filename_check(self, name):
        name = str(name)
        if ' ' in name:
            raise ValueError(f'{name} has spaces in it. \nFiles must not contain spaces')

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
                    self.space_filename_check(target_path)
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
        The elastix image pyramid needs to be specified for each dimension for each resolution.

        This function allows it to specified for just one resolution as it should always be the same for each dimension
        as we are assuming ispotropic data.

        If the the pyramid is already defined per dimesion, it will remain unchanged.


        Example
        -------
        # This pyramid shedule with a single isotropic factor for each stage will be converted to
        MovingImagePyramidSchedule = [6.0, 4.0, 2.0, 1.0, 1.0]

        # Will be converted to
        MovingImagePyramidSchedule = [ 6.0, 6.0, 6.0, 4.0, 4.0, 4.0, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,]

        Parameters
        ----------
        config: dict
            paramters

        """
        pyramid_keys = ['FixedImagePyramidSchedule', 'MovingImagePyramidSchedule']

        for stage in self.config['registration_stage_params']:
            elx_params = stage['elastix_parameters']

            for pk in pyramid_keys:
                if pk in elx_params.keys():

                    if elx_params.get('NumberOfResolutions'):
                        num_res = int(elx_params.get('NumberOfResolutions'))
                    else:
                        raise KeyError('NumberOfResolutions must be specified for each stage')

                    lama_schedule = elx_params[pk]

                    if len(lama_schedule) == num_res:
                        # pyramid is defined one value per stage
                        elastix_shedule = []
                        for i in lama_schedule:
                            elastix_shedule.extend([i, i, i])
                        elx_params[pk] = elastix_shedule

