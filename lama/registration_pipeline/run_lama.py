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

from lama import common
from lama.elastix.invert import InvertLabelMap, InvertMeshes, batch_invert_transform_parameters
from lama.img_processing.normalise import normalise
from lama.img_processing.organ_vol_calculation import label_sizes
from lama.img_processing import glcm3d
from lama.registration_pipeline.validate_config import LamaConfig
from lama.elastix.deformations import make_deformations_at_different_scales
from lama.paths import RegPaths
from lama.qc.metric_charts import make_charts
from lama.elastix.elastix_registration import TargetBasedRegistration, PairwiseBasedRegistration
from lama.staging import staging_metric_maker
from lama.qc.qc_images import make_qc_images_from_config
from lama.stats.standard_stats.data_loaders import DEFAULT_FWHM, DEFAULT_VOXEL_SIZE
from lama.staging.staging_metric_maker import STAGING_METHODS


LOG_FILE = 'LAMA.log'
ELX_PARAM_PREFIX = 'elastix_params_'               # Prefix the generated elastix parameter files
INVERT_CONFIG = 'invert.yaml'
PAD_INFO_FILE = 'pad_info.yaml'


# Set the spacing and origins before registration
SPACING = (1.0, 1.0, 1.0)
ORIGIN = (0.0, 0.0, 0.0)

sys.excepthook = common.excepthook_overide




def run (configfile: Path):
        """
        This is the main function Lama script for generating data from registering volumes
        It reads in the config file, creates directories, and initialises the registration process.

        Looks for paths to inputs relative the directory containing the config file

        Parameters
        ----------
        param config
            A toml config file
        """
        config = LamaConfig(configfile)

        run_registration_schedule(config)




        generate_staging_data(config.staging_method)
        generate_organ_volumes(config)

        create_stats_config: bool = create_stats_config
        config_path: str = configfile

        # The root of the registration project dir
        # self.proj_dir = os.path.dirname(configfile)


        # all paths are relative to the config file directory
        config_dir = os.path.split(os.path.abspath(configfile))[0]

        # logpath = join(self.config_dir, LOG_FILE)
        common.init_logging(config.logpath)

        if __name__ == '__main__':
            logging.info(common.command_line_agrs())

        if not common.test_installation('elastix'):
            raise OSError('Make sure elastix is installed')



        paths = RegPaths(config_dir, config)

        # Number of threads to use during elastix registration
        threads: int = config.get('threads')

        qc: bool
        if not config.get('no_qc'):
            qc_dir = None
        else:
            qc_dir = paths.make('qc')

        outdir = paths.make('output_dir')

        signal.signal(signal.SIGTERM, common.service_shutdown)
        signal.signal(signal.SIGINT, common.service_shutdown)

        memlog_file = Path(outdir) / 'mem_log'
        memmon = common.MonitorMemory(memlog_file)
        memmon.start()

        # The filtype extension to use for registration output, default to nrrd.
        filetype: str = config.filetype

        # Disable QC output
        no_qc: bool = config.get('no_qc')

        logging.info(common.git_log())

        staging_method: str = config.get('staging', 'scaling_factor')

        # Pad the inputs. Also changes the config object to point to these newly padded volumes
        # if config.get('pad_dims'):
        #     self.pad_inputs_and_modify_config()

        logging.info("Registration started")
        final_registration_dir = run_registration_schedule(config)

        if config.get('glcm'):
            create_glcms()

        if config.get('skip_transform_inversion'):
            logging.info('Skipping inversion of transforms')
        else:
            make_inversion_transform_files()

            invert_status = invert_volumes(config)

            if invert_status:

                generate_organ_volumes(config)

        generate_staging_data(staging_method)

        if not no_qc:
            registered_midslice_dir = Path(paths.make('registered_midslice_dir', parent=qc_dir))
            inverted_label_overlay_dir = Path(paths.make('inverted_label_overlay_dir', parent=qc_dir))
            make_qc_images_from_config(config, Path(outdir), Path(registered_midslice_dir), Path(inverted_label_overlay_dir))

        memmon.stop()
        memmon.join()


    #     self.plot_memory_usage(memlog_file)
    #
    # def plot_memory_usage(self, memory_file):
    #     df = pd.read_csv(memory_file)
    #     df.time = df.time.astype(pd.datetime)
    #     sns.lineplot(x='time', y='mem %', data=df)
    #     plt.savefig(Path(self.outdir) / 'mem_log.png')


def generate_staging_data(config, staging_method):
    """
    Generate staging data from the registration resutls

    Current options are
        scaling_factor (default)
        label_length
        whicle_volme

    Parameters
    ----------
    staging_method: str
        the staging method to use
    """

    if staging_method == 'scaling_factor':
        logging.info('Doing stage estimation - scaling factor')
        stage_dir = config.get_affine_or_similarity_stage_dir()
        staging_metric_maker.scaling_factor_staging(stage_dir, self.outdir)

    elif staging_method == 'embryo_volume':
        logging.info('Doing stage estimation - whole embryo volume')
        # Get the first registration stage dir
        first_stage_reg_dir = config['registration_stage_params'][0]['stage_id']
        inv_mask_path = config['inverted_stats_masks']
        staging_metric_maker.whole_volume_staging(inv_mask_path, self.outdir)

def get_affine_or_similarity_stage_dir(self):
    """
    Get the output path to the first occurence of a similarity or affine registration
    Returns
    -------
    str: path to siilarity or affine registration dir
    """
    for stage_info in self.config['registration_stage_params']:
        if stage_info['elastix_parameters']['Transform'] in ('SimilarityTransform', 'AffineTransform'):
            reg_dir = join(self.outdir, 'registrations', stage_info['stage_id'])
            return reg_dir

def make_inversion_transform_files(self):
    """
    Create inversion transform parameter files that can be used to invert volumes in population average space back
    onto the inputs


    Parameters
    ----------
    config

    Returns
    -------

    """
    logging.info('inverting transforms')
    tform_invert_dir = self.paths.make('inverted_transforms')

    # Path to create a config that specifies the orrder of inversions
    self.invert_config = join(tform_invert_dir, INVERT_CONFIG)

    batch_invert_transform_parameters(self.config_path, self.invert_config, tform_invert_dir, self.threads)

def invert_volumes(self, config):
    """
    Invert volumes, such as masks and labelmaps from population average space to input volumes space using
    pre-calculated elastix inverse transform parameter files
    Parameters
    ----------
    config

    Returns
    -------

    """

    # 240918
    # These tweo lines are duplicated from make_inversion_tranform_file. Just a bodge for now so can use restart_at_stage == invert_volumes
    tform_invert_dir = self.paths.make('inverted_transforms')
    self.invert_config = join(tform_invert_dir, INVERT_CONFIG)


    status = True
    if config.get('stats_mask'):
        mask_path = self.paths['stats_mask']
        self.invert_labelmap(mask_path, 'inverted_stats_masks')
    else:
        status = False

    if config.get('label_map'):
        labelmap = self.paths['label_map']
        self.invert_labelmap(labelmap, 'inverted_labels')
    else:
        status = False

    return status

def generate_organ_volumes(self):

    # Get the final inversion stage
    with open(self.invert_config, 'r') as fh:
        first_stage = yaml.load(fh)['inversion_order'][-1]

    inverted_label_dir =  join(self.paths.get('inverted_labels'), first_stage)
    # inverted_mask_dir = join(self.paths.get('inverted_stats_masks'), first_stage)  030918 not using normalised volumes
    out_path = self.paths.get('organ_vol_result_csv')
    label_sizes(inverted_label_dir, out_path)

def invert_labelmap(self, label_file: Path, name=None):
    """
    invert a labelmap.
    Parameters
    ----------
    label_file:
        labelmap or mask file
    name:
        name to call inversion directory

    """

    if not os.path.isfile(label_file):
        logging.info('labelmap: {} not found')
        return

    label_inversion_dir = self.paths.make(name)

    ilm = InvertLabelMap(self.invert_config, label_file, label_inversion_dir, threads=self.threads)
    ilm.run()
    return label_inversion_dir

def invert_isosurfaces(self):
    """
    Invert a bunch of isosurfaces that were proviously generated from the target labelmap
    For this to work, the target needs to be the largest of all volumes used in registration, as padding the target
    will through the mesh and taget corrdinates out of sync
    :return:
    """
    if not self.config.get('isosurface_dir'):
        logging.info('isosurface directory not in config file')
        return

    # TODO: put a check for target being the largest volume ?

    mesh_dir = self.paths['isosurface_dir']
    if not os.path.isdir(mesh_dir):
        logging.info('Mesh directory: {} not found'.format(mesh_dir))

    iso_out = self.paths.make('inverted_isosurfaces')

    logging.info('mesh inversion started')
    for mesh_path in common.get_file_paths(mesh_dir):
        im = InvertMeshes(self.invert_config, mesh_path, iso_out)
        im.run()

def run_registration_schedule(config):
    """
    Run the registrations specified in the config file
    :param config: Config dictionary
    :return: 0 for success, or an error message
    """

    stage_params = generate_elx_parameters(config, config['pairwise_registration'])

    # Make dir to put averages in
    avg_dir = self.paths.make('averages')

    root_reg_dir = self.paths.make('root_reg_dir', mkdir='force')

    # if True: create a new fixed volume by averaging outputs
    # if False: use the same fixed volume at each registration stage
    regenerate_target = config.get('generate_new_target_each_stage')
    if regenerate_target:
        logging.info('Creating new target each stage for population average creation')
    else:
        logging.info('Using same target for each stage')

    filetype = config['filetype']

    # Set the moving vol dir and the fixed image for the first stage
    moving_vols_dir = config['inputs']

    if not self.no_qc:
        input_histogram_dir = self.paths.make('input_image_histograms', parent=self.qc_dir)
        make_histograms(moving_vols_dir, input_histogram_dir)

    fixed_vol = self.paths['fixed_volume']

    # Create a folder to store mid section coronal images to keep an eye on registration process

    if not self.no_qc:
        qc_image_dir = self.paths.make('qc_registered_images', parent=self.qc_dir)
        qc_metric_dir = self.paths.make('metric_charts', parent=self.qc_dir)

    reg_stages = config['registration_stage_params']

    for i, reg_stage in enumerate(reg_stages):

        tform_type = reg_stage['elastix_parameters']['Transform']
        euler_stage = True if tform_type == 'EulerTransform' else False
        # affine_similarity_stage = True if tform_type in ['AffineTransform', 'SimilarityTransform'] else False

        if do_pairwise and not euler_stage:
            logging.info('doing pairwise registration')
            RegMethod = PairwiseBasedRegistration

        elif not do_pairwise:
            logging.info('using target-based registration')
            RegMethod = TargetBasedRegistration

        elif do_pairwise and euler_stage:
            RegMethod = TargetBasedRegistration
            logging.info(
                'using target-based registration for initial rigid stage of pairwise registrations')

        #  Make the stage output dir
        stage_id = reg_stage['stage_id']
        stage_dir = join(root_reg_dir, stage_id)

        common.mkdir_force(stage_dir)

        logging.info("### Current registration step: {} ###".format(stage_id))

        # Make the elastix parameter file for this stage
        elxparam = stage_params[stage_id]
        elxparam_path = join(stage_dir, ELX_PARAM_PREFIX + stage_id + '.txt')

        with open(elxparam_path, 'w') as fh:
            if elxparam:  # Not sure why I put this here
                fh.write(elxparam)

        # TODO: no fixed mask after 1st stage if doing pairwise
        # For the first stage, we can use the fixed mask for registration.
        # Sometimes helps with the 'too many samples map outside fixed image' problem
        if i < 2:  # TODO: shall we keep the fixed mask throughout?
            fixed_mask = self.paths['fixed_mask']

        # # If we are doing target-based phenotype detection, we can used the fixed mask for every stage
        # elif not self.config.get('generate_new_target_each_stage'):
        #     fixed_mask = self.config.get('fixed_mask')
        else:
            fixed_mask = None

        # if fixed_mask:
        #     fixed_mask = join(self.config_dir, fixed_mask)


        # Do the registrations
        registrator = RegMethod(elxparam_path,
                                moving_vols_dir,
                                stage_dir,
                                self.paths,
                                self.filetype,
                                self.threads,
                                fixed_mask,
                                self.config_dir
                                )
        if (not do_pairwise) or (do_pairwise and euler_stage):
            registrator.set_target(fixed_vol)

        registrator.run()

        # Make average
        average_path = join(avg_dir, '{0}.{1}'.format(stage_id, filetype))
        registrator.make_average(average_path)

        if not self.no_qc:

            stage_metrics_dir = join(qc_metric_dir, stage_id)
            self.paths.make(stage_metrics_dir)
            make_charts(stage_dir, stage_metrics_dir)

        # Reregister the inputs to this stage to the average from this stage to create a new average
        if config.get('re-register_each_stage'):
            average_path, stage_dir = self.repeat_registration(average_path, moving_vols_dir, elxparam_path, average_path)

        # Setup the fixed, moving and mask_dir for the next stage, if there is one
        if i + 1 < len(config['registration_stage_params']):
            if regenerate_target:
                fixed_vol = average_path  # The avergae from the previous step

            moving_vols_dir = stage_dir  # The output dir of the previous registration

    if config.get('generate_deformation_fields'):
        make_vectors = not config.get('skip_deformation_fields')
        make_deformations_at_different_scales(config, root_reg_dir, self.outdir, make_vectors, self.threads,
                                              filetype=config.get('filetype'), skip_histograms=config.get('no_qc'))

    logging.info("### Registration finished ###")
    return stage_dir  # Return the path to the final registrerd images

    return target, moving_vols_dir


def create_glcms(self):
    """
    Create grey level co-occurence matrices. This is done in the main registration pipeline as we don't
    want to have to create GLCMs for the wildtypes multiple times when doing phenotype detection
    """
    logging.info("Creating GLCMS")
    glcm_out_dir = self.paths.make('glcms')  # The vols to create glcms from
    if self.config.get('fixed_mask'):
        mask_path = self.paths['fixed_mask']
        mask = common.img_path_to_array(mask_path)
    else:
        logging.warn("Cannot make GLCMs without a mask")
        return
    glcm3d.pyradiomics_glcm(self.final_registration_dir, glcm_out_dir, mask)

    logging.info("Finished creating GLCMs")

def normalise_registered_images(self, stage_dir, norm_dir, norm_roi):

    roi_starts = norm_roi[0]
    roi_ends = norm_roi[1]


    logging.info('Normalised registered output using ROI: {}:{}'.format(
        ','.join([str(x) for x in roi_starts]), ','.join([str(x) for x in roi_ends])))
    normalise(stage_dir, norm_dir, roi_starts, roi_ends)

def generate_elx_parameters(self, config, do_pairwise=False):
    """
    Generate a dict of elestix parameters for each stage. Merge global paramters into each stage. inherit
    parameters from previous stages when required, and overide with stage-specific parameters

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
    """
    stage_params = OrderedDict()
    for i, reg_stage in enumerate(config['registration_stage_params']):

        stage_id = reg_stage['stage_id']
        elxparams_formated = []
        param_vals = []  # The parameters to format

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
        param_vals.extend([global_params, reg_stage['elastix_parameters']])

        for p in param_vals:
            for param_name, param_value in p.items():
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

