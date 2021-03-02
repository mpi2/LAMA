#! /usr/bin/env python3

"""
=====================================================
MRC Harwel Registration pipeline
=====================================================
..module:: The Medical Research Council registration pipeline
:synopsis: Pipeline to create embryo average volumes.
.. moduleauthor:: Neil Horner


This is the main module of the registration of the mouse embryo registration pipeline.
It can be used for creating population averages and for the creation of data to input into the lama stats pipeline

Usage
-----
An example config file can be found :download:`here<example_configs/pop_average_config.yaml>`.

Config file
-----------

* The config file is in yaml format (http://yaml.org/)
* Paths are relative to the directory that contains the config file
* Everything after a '#' a comment
* The following entries should appear at the top of the config file:


.. code-block:: yaml

    fixed_volume: target/2608145_deformable_avg_8bit_28um.nrrd # The fixed volume
    inputs: inputs  # directory with moving volumes
    fixed_mask: target/fixed_mask_for_reg.nrrd  # fixed mask path
    pad_dims: true # Pads all the volumes so all are the same dimensions. Finds the largest dimension from each volume
    pad_dims: [300, 255, 225]  # this specifies the dimensions tyo pad to
    threads: 10  # number of cpu cores to use
    filetype: nrrd  # the output file format - nrrd, nii, tiff
    compress_averages: true  # compress the averages to save disk space
    generate_new_target_each_stage: true  # true for creating an average. false for phenotype detection
    
    staging entry. this allows for the automatoc determination of stage using various surrogates
    staging: scaling_factor
  
    staging: whole_volume  # not yet implemented
    staging_volume: path/to/mask

    staging: label_length
    staging_volume: path/to/label


The next section specifies parameters to be used by elastix at all stages of registration


.. code-block:: yaml


    global_elastix_params:
       FixedInternalImagePixelType: short
       MovingInternalImagePixelType: short
       FixedImageDimension: 3

Then everything under the following section species each stage of the multi-resolution registration

.. code-block:: yaml

    registration_stage_params:

Then we specify a registration stage like this

.. code-block:: yaml

    - stage_id: rigid
      elastix_parameters:
        Metric: AdvancedNormalizedCorrelation
        Registration: MultiResolutionRegistration

| If generate_new_target_each_stage: true
| Then the average from this stage will be used in the next stage
|
| We can add some other parameters before the elastix parameters to the stage settings as below


.. code-block:: yaml

    - stage_id: deformable_to_8
      inherit_elx_params: deformable_to_128
      do_analysis: true
      normalise_registered_output: [[200, 240, 230], [210, 250, 240]]

inherit_elx_params: This takes the elastix paramteters from the named stage. Any elastix parameters specified
after this will overide the inherited parameters

do_analysis: if set to true will generate deformation fields and spatial jacobians from this registration stage

normalise_registered_output: in order to do intensity-based analsysis on the registered output, it should be
normalised first. The format [[start indices, [end indices]] specifies an ROI from an area of background in the
target. This average of this region in the outputs will be used as as the new zero value for the image


"""
from lama import common
common.disable_warnings_in_docker()

import os
import copy
from logzero import logger as logging
from collections import OrderedDict
from os.path import join
import SimpleITK as sitk
import yaml
import sys
from pathlib import Path
import signal
import shutil

from lama.elastix.invert_volumes import InvertLabelMap, InvertMeshes
from lama.elastix.invert_transforms import batch_invert_transform_parameters
from lama.elastix.reverse_registration import reverse_registration
from lama.img_processing.organ_vol_calculation import label_sizes
from lama.img_processing import glcm3d
from lama.registration_pipeline.validate_config import LamaConfig, LamaConfigError
from lama.elastix.deformations import make_deformations_at_different_scales
from lama.qc.metric_charts import make_charts
from lama.elastix.elastix_registration import TargetBasedRegistration, PairwiseBasedRegistration
from lama.staging import staging_metric_maker
from lama.qc.qc_images import make_qc_images
from lama.qc.folding import folding_report
from lama.stats.standard_stats.data_loaders import DEFAULT_FWHM, DEFAULT_VOXEL_SIZE
from lama.elastix import INVERT_CONFIG, REG_DIR_ORDER_CFG
from lama.monitor_memory import MonitorMemory
from lama.common import cfg_load
from lama.segmentation_plugins import plugin_interface

LOG_FILE = 'LAMA.log'
ELX_PARAM_PREFIX = 'elastix_params_'               # Prefix the generated elastix parameter files
PAD_INFO_FILE = 'pad_info.yaml'

temp_debug_fix_folding = True


# Set the spacing and origins before registration
SPACING = (1.0, 1.0, 1.0)
ORIGIN = (0.0, 0.0, 0.0)

sys.excepthook = common.excepthook_overide


def run(configfile: Path):
        """
        This is the main function Lama script for generating data from registering volumes
        It reads in the config file, creates directories, and initialises the registration process.

        Looks for paths to inputs relative the directory containing the config file

        Parameters
        ----------
        param config
            A toml config file
        """
        try:
            config = LamaConfig(configfile)
        except OSError as e:
            logging.error(f'Cannot open LAMA config file: {str(configfile)}\n{e}')
            raise
        except Exception as e:
            raise(LamaConfigError(e))

        config.mkdir('output_dir')
        qc_dir = config.mkdir('qc_dir')
        config.mkdir('average_folder')
        config.mkdir('root_reg_dir')

        # TODO find the histogram batch code
        # if not config['no_qc']:
        #     input_histogram_dir = config.mkdir('input_image_histograms')
        #     make_histograms(config['inputs'], input_histogram_dir)

        logpath = config.config_path.parent / LOG_FILE  # Make log in same directory as config file
        common.init_logging(logpath)

        if not common.test_installation('elastix'):
            raise OSError('Make sure elastix is installed')

        # Catch ctr-c signals so we can write that to logs
        # signal.signal(signal.SIGTERM, common.service_shutdown)
        signal.signal(signal.SIGINT, common.service_shutdown)

        mem_monitor = MonitorMemory(Path(config['output_dir']).absolute())

        # Disable QC output?
        no_qc: bool = config['no_qc']

        logging.info(common.git_log())  # If running from a git repo, log the branch and commit #

        logging.info("Registration started")

        first_stage_only = config['skip_forward_registration']
        # If we only want the reverse label propagation we just need the initial rigid registration to act as the
        # Fixed image for the moving populaiton average
        final_registration_dir = run_registration_schedule(config, first_stage_only=first_stage_only)

        if not first_stage_only:
            neg_jac = make_deformations_at_different_scales(config)
            folding_report(neg_jac, config['output_dir'], config['label_info'], outdir=config['output_dir'])

            create_glcms(config, final_registration_dir)

        # Write out the names of the registration dirs in the order they were run
        with open(config['root_reg_dir'] / REG_DIR_ORDER_CFG, 'w') as fh:
            for reg_stage in config['registration_stage_params']:
                fh.write(f'{reg_stage["stage_id"]}\n')
                if first_stage_only:
                    break

        if config['skip_transform_inversion']:
            logging.info('Skipping inversion of transforms')
        else:
            logging.info('inverting transforms')

            if config['label_propagation'] == 'reverse_registration':
                reverse_registration(config)
            else:  # invert_transform method is the default
                batch_invert_transform_parameters(config)

            logging.info('inverting volumes')
            invert_volumes(config)

            # Now that labels have been inverted, should we delete the transorm files?
            if config['delete_inverted_transforms']:
                shutil.rmtree(config['output_dir'] / 'inverted_transforms')

            if config['label_map']:

                generate_organ_volumes(config)

                if config['seg_plugin_dir']:
                    plugin_interface.secondary_segmentation(config)

        if not generate_staging_data(config):
            logging.warning('No staging data generated')

        if not no_qc:

            rev_reg = True if config['label_propagation'] == 'reverse_registration' else False
            make_qc_images(config.config_dir, config['fixed_volume'], qc_dir, mask=None,
                           reverse_reg_propagation=rev_reg)

        mem_monitor.stop()

        return True

# def get_whole_embryo_mask(config: LamaConfig):
#     mask_root = config['out_dir'] / 'inverted_stats_masks'
#     if config['label_propagation'] == 'reverse_registration':
#         mask ='c'
#     elif config['label_propagation'] == 'invert_transform':
#         mask = config['inverted_stats_masks'] / 'rigid'

def generate_staging_data(config: LamaConfig):
    """
    Generate staging data from the registration results
    """

    staging_method = config['staging']

    if staging_method == 'scaling_factor':
        logging.warn('Scaling factor not working at the moment')
        return
        logging.info('Doing stage estimation - scaling factor')
        stage_dir = get_affine_or_similarity_stage_dir(config)
        staging_metric_maker.scaling_factor_staging(stage_dir, config['output_dir'])
        return True

    elif staging_method == 'embryo_volume':
        logging.info('Doing stage estimation - whole embryo volume')

        # Get the dir name of the stage that we want to calculate organ volumes from (rigid, affine)
        stage_to_get_volumes = get_affine_or_similarity_stage_dir(config)
        if not stage_to_get_volumes:
            logging.warn('Cannot find a similarity or affine stage to generate whole embryo volume staging data.')
            return

        logging.info('Generating whole embryo volume staging data')
        # Get the root of the inverted masks for thie current specimen
        inv_mask_root = config['inverted_stats_masks']

        inv_mask_stage_dir = inv_mask_root / stage_to_get_volumes.name
        staging_metric_maker.whole_volume_staging(inv_mask_stage_dir, config['output_dir'])
        return True


def get_affine_or_similarity_stage_dir(config: LamaConfig) -> Path:
    """
    Get the output path to either a similarity or affine registration.
    Earliest stage in pipeline takes precedent if multiple valid stages present.

    Returns
    -------
    path to similarity or affine registration dir
    """

    for stage_info in config['registration_stage_params']:
        if stage_info['elastix_parameters']['Transform'] in ('SimilarityTransform', 'AffineTransform'):
            reg_dir = config['root_reg_dir'] / stage_info['stage_id']

            return reg_dir


def invert_volumes(config: LamaConfig):
    """
    Invert volumes, such as masks and labelmaps from population average space to input volumes space using
    pre-calculated elastix inverse transform parameter files

    Returns
    -------
    Status of inverions for masks and labels

    """

    invert_config = config['inverted_transforms'] / INVERT_CONFIG

    if config['stats_mask']:
        mask_inversion_dir = config.mkdir('inverted_stats_masks')
        InvertLabelMap(invert_config, config['stats_mask'], mask_inversion_dir, threads=config['threads']).run()

    if config['label_map']:
        labels_inverion_dir = config.mkdir('inverted_labels')
        InvertLabelMap(invert_config, config['label_map'], labels_inverion_dir, threads=config['threads']).run()


def generate_organ_volumes(config: LamaConfig):

    # Get the final inversion stage
    invert_config = config['inverted_transforms'] / INVERT_CONFIG

    first_stage = cfg_load(invert_config)['inversion_order'][-1]

    inverted_label_dir =  config['inverted_labels'] / first_stage

    out_path = config['organ_vol_result_csv']

    # Generate the organ volume csv
    label_sizes(inverted_label_dir, out_path)


# def invert_isosurfaces(self):
#     """
#     Invert a bunch of isosurfaces that were proviously generated from the target labelmap
#     For this to work, the target needs to be the largest of all volumes used in registration, as padding the target
#     will through the mesh and taget corrdinates out of sync
#     :return:
#     """
#     if not self.config.get('isosurface_dir'):
#         logging.info('isosurface directory not in config file')
#         return
#
#     # TODO: put a check for target being the largest volume ?
#
#     mesh_dir = self.paths['isosurface_dir']
#     if not os.path.isdir(mesh_dir):
#         logging.info('Mesh directory: {} not found'.format(mesh_dir))
#
#     iso_out = self.paths.make('inverted_isosurfaces')
#
#     logging.info('mesh inversion started')
#     for mesh_path in common.get_file_paths(mesh_dir):
#         im = InvertMeshes(self.invert_config, mesh_path, iso_out)
#         im.run()


def run_registration_schedule(config: LamaConfig, first_stage_only=False) -> Path:
    """
    Run the registrations specified in the config file

    Parameters
    ----------
    config: Parsed and validated lama config
    first_stage_only: If True, just do the initial rigid stage

    Returns
    -------
    The path to the final registrered images
    """
    st = config['stage_targets']
    if st:
        with open(st, 'r') as stfh:
            stage_targets = cfg_load(stfh)['targets']
        if len(config['registration_stage_params']) != len(stage_targets):
            logging.error(f'Len stage targets: {len(stage_targets)}')
            logging.error(f'Len reg stages: {len(config["registration_stage_params"])}')
            raise LamaConfigError("restage len != number of registration stages")

    # Create a folder to store mid section coronal images to keep an eye on registration process
    if not config['no_qc']:
        qc_metric_dir = config['metric_charts_dir']

    elastix_stage_parameters = generate_elx_parameters(config, do_pairwise=config['pairwise_registration'])
    regenerate_target = config['generate_new_target_each_stage']

    if regenerate_target and st:
        raise LamaConfigError('cannot have regenerate_target and stage_targets')

    if regenerate_target:
        logging.info('Creating new target each stage for population average creation')
    else:
        logging.info('Using same target for each stage')

    # Set the moving volume dir and the fixed image for the first stage
    moving_vols_dir = config['inputs']

    # Set the fixed volume up for the first stage. This will checnge each stage if doing population average
    if st:
        fixed_vol = stage_targets[0]
    else:
        fixed_vol = config['fixed_volume']

    for i, reg_stage in enumerate(config['registration_stage_params']):

        tform_type = reg_stage['elastix_parameters']['Transform']
        euler_stage = True if tform_type == 'EulerTransform' else False
        # affine_similarity_stage = True if tform_type in ['AffineTransform', 'SimilarityTransform'] else False

        if config['pairwise_registration']:
            if not euler_stage:
                logging.info('doing pairwise registration')
                reg_method = PairwiseBasedRegistration
            else:
                reg_method = TargetBasedRegistration
                logging.info(
                    'using target-based registration for initial rigid stage of pairwise registrations')

        else:
            logging.info('using target-based registration')
            reg_method = TargetBasedRegistration

        #  Make the stage output dir
        stage_id = reg_stage['stage_id']
        stage_dir = config.stage_dirs[stage_id]

        common.mkdir_force(stage_dir)

        logging.info("### Current registration step: {} ###".format(stage_id))

        # Make the elastix parameter file for this stage
        elxparam = elastix_stage_parameters[stage_id]
        elxparam_path = join(stage_dir, ELX_PARAM_PREFIX + stage_id + '.txt')

        with open(elxparam_path, 'w') as fh:
            if elxparam:  # Not sure why I put this here
                fh.write(elxparam)

        # if i < 2:  # TODO: shall we keep the fixed mask throughout? I think we should in next release
        #     fixed_mask = config['fixed_mask']

        # # If we are doing target-based phenotype detection, we can used the fixed mask for every stage
        if not config['generate_new_target_each_stage']:
            fixed_mask = config['fixed_mask']
        else:
            fixed_mask = None

        # Do the registrations
        registrator = reg_method(elxparam_path,
                                 moving_vols_dir,
                                 stage_dir,
                                 config['filetype'],
                                 config['threads'],
                                 fixed_mask
                                 )

        if (not config['pairwise_registration']) or (config['pairwise_registration'] and euler_stage):
            registrator.set_target(fixed_vol)

        if reg_stage['elastix_parameters']['Transform'] == 'BSplineTransform':
            logging.info(f'Folding correction for stage {stage_id} set')
            registrator.fix_folding = config['fix_folding']  # Curently only works for TargetBasedRegistration

        registrator.run()  # Do the registrations for a single stage

        # Make average from the stage outputs
        if regenerate_target:
            average_path = join(config['average_folder'], '{0}.{1}'.format(stage_id, config['filetype']))
            registrator.make_average(average_path)

        if not config['no_qc']:

            stage_metrics_dir = qc_metric_dir / stage_id
            common.mkdir_force(stage_metrics_dir)
            make_charts(stage_dir, stage_metrics_dir)

        # Setup the fixed and moving for the next stage, if there is one
        if i + 1 < len(config['registration_stage_params']):
            if regenerate_target:
                fixed_vol = average_path  # The avergae from the previous step
            elif st:
                fixed_vol = stage_targets[i+1]

            moving_vols_dir = stage_dir  # Set the output of the current stage top be the input of the next

        if first_stage_only:
            return stage_dir

    logging.info("### Registration finished ###")

    return stage_dir


def create_glcms(config: LamaConfig, final_reg_dir):
    """
    Create grey level co-occurence matrices. This is done in the main registration pipeline as we don't
    want to have to create GLCMs for the wildtypes multiple times when doing phenotype detection.
    """
    if not config['glcm']:
        return
    logging.info("Creating GLCMS")
    glcm_dir = config.mkdir('glcm_dir')

    mask_path = config['fixed_mask']
    if mask_path:
        mask = common.img_path_to_array(mask_path)
    else:
        logging.warn("Cannot make GLCMs without a mask")
        return

    glcm3d.pyradiomics_glcm(final_reg_dir, glcm_dir, mask)
    logging.info("Finished creating GLCMs")


def generate_elx_parameters(config: LamaConfig,
                            do_pairwise: bool = False,
                            ) -> OrderedDict:
    """
    Generate an ordered dictionary of elastix parameters for each stage.
    Merge global parameters into each stage.
    Inherit parameters from previous stages when required (if inherit_elx_params is specified)
    Overide globals parameters with stage-specific parameters.

    Parameters
    ----------
    config: dict
        the main config
    do_pairwise: Bool
        if True set elestix parameters to write result image

    Returns
    -------
    dict:
        elastix parameters for each stage
        stage_id: {elastix paramters}...
    """
    stage_params = OrderedDict()
    for i, reg_stage in enumerate(config['registration_stage_params']):

        stage_id = reg_stage['stage_id']
        elxparams_formated = []
        # param_vals = []  # The parameters to format
        params = {}


        inherit_stage = reg_stage.get('inherit_elx_params')
        if inherit_stage:  # Use another stage's params as a starting point

            # Find the referenced stage's parameters
            referenced_stage_found = False
            for stage in config['registration_stage_params']:
                if stage['stage_id'] == inherit_stage:
                    referenced_stage_found = True
                    inherited_params = copy.deepcopy(stage['elastix_parameters'])
                    #  Overwrite the inherited params with this stage-specific ones
                    for k, v in reg_stage['elastix_parameters'].items():
                        inherited_params[k] = v
                    # Now replace the stage-specific params with the merged new param dict
                    reg_stage['elastix_parameters'] = inherited_params

            if not referenced_stage_found:
                sys.exit("Cannot inherit parameters from the stage: {}. Stage Not found. Check the config file".format(
                    stage_params['inherit_elx_params']))

        # Merge the global and stage-specific parameters
        global_params = config['global_elastix_params']

        # Make sure result image is written if using taregt based registration. Needed for the moving images of the
        # next stage
        if (do_pairwise and reg_stage['elastix_parameters']['Transform'] == 'EulerTransform') or (not do_pairwise):
            reg_stage['elastix_parameters']['WriteResultImage'] = 'true'
        else:
            reg_stage['elastix_parameters']['WriteResultImage'] = 'false'
        global_params.pop('WriteResultImage', None)

        # global_params.update(extra_global_params)
        params.update(reg_stage['elastix_parameters'])
        params.update(global_params)

        # for p in param_vals:
        for param_name, param_value in params.items():
            if isinstance(param_value, list):  # parameter with multiple values for each resolution
                # check whether we have didgets or strings, the latter need quoting
                if all(common.is_number(v) for v in param_value):
                    val_list = [str(s).format() for s in param_value]
                    val_string = ' '.join(val_list)
                else:
                    val_list = [str(s).format() for s in param_value]
                    val_quoted = ['"{}"'.format(s) for s in val_list]
                    val_string = ' '.join(val_quoted)
                elxparams_formated.append('({0}  {1})\n'.format(param_name, val_string))
            else:  # Find a nicer way to do this
                if common.is_number(param_value):
                    elxparams_formated.append('({0}  {1})\n'.format(param_name, param_value))
                else:
                    elxparams_formated.append('({0}  "{1}")\n'.format(param_name, param_value))

        stage_params[stage_id] = ''.join(elxparams_formated)
    return stage_params


def write_stats_config(config_path: Path, config, key_values):
    """
    Write a config for input into the stats pipeline
    """
    out_path = config_path.parent / 'stats.yaml'
    keys = list(key_values.keys())

    # TODo make the name in the values not paths
    stats_cfg = dict(
        # Put in target dir
        stats_types = ['intensity', 'jacobians', 'organ_volume'],
        shape=config['pad_dims'],
        mask=config['stats_mask'],
        label_map=config['label_map'],
        blur_fwhm=DEFAULT_FWHM,
        voxel_size=DEFAULT_VOXEL_SIZE
    )

    with open(out_path, 'w') as fh:
        fh.write(yaml.dump(stats_cfg))


def set_origins_and_spacing(volpaths):
    for vol in volpaths:
        im = sitk.ReadImage(vol)
        if im.GetSpacing() != SPACING or im.GetOrigin() != ORIGIN:
            im.SetSpacing(SPACING)
            im.SetOrigin(ORIGIN)
            sitk.WriteImage(im, vol)


def make_histograms(in_dir, out_dir):
    """
    Make histograms of a series of input images and output to a html file
    """
    return
    try:
        hist_batch(in_dir, out_dir, remove_zeros=True)
    except Exception as e:
        logging.exception('There was an error generating the histograms\n{}'.format(e))


def mkdir_if_not_exists(dir_):
    if not os.path.exists(dir_):
        os.makedirs(dir_)

