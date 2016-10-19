#!/usr/bin/python

"""
=====================================================
MRC Harwel Registration pipeline
=====================================================
..module:: The Medical Research Council registration pipeline
:synopsis: Pipeline to create embryo average volumes.
.. moduleauthor:: Neil Horner


This is the main module of the registration of the mouse embryo registration pipeline.
It can be used for creating population averages and for the detection of anatomical phenotypes ( see pheno_detect.py)


Requirements
------------
* Python 2.7
* elastix (http://elastix.isi.uu.nl/) for doing the registrations. Works with v4.7
* PyYAML
* SimpleITK

Usage
-----

.. code-block:: bash

    ./mrch_regpipeline.py -c wildtype_config.yaml


An example config file can be found :download:`here<example_configs/pop_average_config.yaml>`.

Config file
-----------


* The config file is in yaml format (http://yaml.org/)
* Paths are relative to the directory that contains the config file
* Everything after a '#' a comment
* The following entries should appear at the top of the config file:


.. code-block:: yaml

    output_dir: out  # Registration output goes here
    fixed_volume: target/2608145_deformable_avg_8bit_28um.nrrd # The fixed volume
    inputs: inputs  # directory with moving volumes
    fixed_mask: target/fixed_mask_for_reg.nrrd  # fixed mask path
    pad_dims: true # Pads all the volumes so all are the same dimensions. Finds the largest dimension from each volume
    pad_dims: [300, 255, 225]  # this specifies the dimensions tyo pad to
    threads: 10  # number of cpu cores to use
    filetype: nrrd  # the output file format - nrrd, nii, tiff
    compress_averages: true  # compress the averages to save disk space
    generate_new_target_each_stage: true  # true for creating an average. false for phenotype detection

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

import os
from os.path import join, splitext, basename, relpath
import subprocess
import shutil
import sys
import argparse
from collections import OrderedDict, defaultdict
import copy
import itertools
import logging
import yaml
import SimpleITK as sitk
import numpy as np
import scipy.misc
from normalise import normalise
from invert import InvertLabelMap, InvertMeshes, batch_invert_transform_parameters
import common
try:
    import glcm3d
except ImportError:
    glcm3d = False
from calculate_organ_size import calculate_volumes
from validate_config import validate_reg_config
from deformations import make_deformations_at_different_scales
from paths import RegPaths
from metric_charts import make_charts
from elastix_registration import TargetBasedRegistration, PairwiseBasedRegistration
from utilities.histogram_batch import batch as hist_batch
from roi_overlay import make_normalization_roi_qc_images
from pad import pad_volumes

LOG_FILE = 'LAMA.log'
ELX_PARAM_PREFIX = 'elastix_params_'               # Prefix the generated elastix parameter files
INDV_REG_METADATA = 'reg_metadata.yaml'            # file name  for singleregistration metadata file
INVERT_CONFIG = 'invert.yaml'
ORGAN_VOLS_OUT = 'organ_volumes.csv'
PAD_INFO_FILE = 'pad_info.yaml'


# Set the spacing and origins before registration
SPACING = (1.0, 1.0, 1.0)
ORIGIN = (0.0, 0.0, 0.0)

SINGLE_THREAD_METRICS = ['TransformRigidityPenalty']


class RegistraionPipeline(object):
    def __init__(self, configfile, create_modified_config=True):
        """This is the main function that is called by the GUI or from the command line.
        Reads in the config file, Creates directories, and initialises the registration process

        :param config: either a yaml file from the command line. Or a dict from the GUI
        :param: callback, function to return messages back to the GUI
        :return: nothing yet
        """

        # Log all uncaught exceptions
        sys.excepthook = common.excepthook_overide

        self.create_modified_config = create_modified_config
        self.config_path = configfile

        # The root of the registration project dir
        self.proj_dir = os.path.dirname(configfile)

        # Get the contents of the config file into a dict
        self.config = config = parse_yaml_config(configfile)

        # all paths are relative to the config file directory
        self.config_dir = os.path.split(os.path.abspath(configfile))[0]

        logpath = join(self.config_dir, LOG_FILE)
        common.init_logging(logpath)

        self.test_elastix_installation()

        # Validate the config file to look for common errors. Add defaults
        validate_reg_config(config, self.config_dir)

        self.paths = RegPaths(self.config_dir, self.config)

        self.outdir = self.paths.make('outdir')

        # Number oif threads to use during elastix registration
        threads = self.config.get('threads')
        if threads:
            self.threads = str(threads)

        # The filtype extension to use for registration output, use nrrd if not set
        self.filetype = config.get('filetype', 'nrrd')

        # Disable QC output
        self.no_qc = config.get('no_qc')

        logging.info(common.git_log())

        self.restart_stage = config.get('restart_at_stage')

        # Pad the inputs. Also changes the config object to point to these newly padded volumes

        if config.get('pad_dims') and not self.restart_stage:
            self.pad_inputs_and_modify_config()

        logging.info("Registration started")
        self.run_registration_schedule(config)

        if self.config.get('skip_transform_inversion'):
            logging.info('Skipping inversion of transforms')
        else:
            logging.info('inverting transforms')
            tform_invert_dir = self.paths.make('inverted_transforms')

            self.invert_config = join(tform_invert_dir, INVERT_CONFIG)
            self._invert_elx_transform_parameters(tform_invert_dir)

            if config.get('label_map_path'):
                self.invert_labelmap()

            if self.config.get('isosurface_dir'):
                self.invert_isosurfaces()

    def _invert_elx_transform_parameters(self, invert_out):
        """
        Invert the elastix output transform parameters. The inverted parameter files can then be used for inverting
        labelmaps and statistics overlays etc.
        """

        batch_invert_transform_parameters(self.config_path, self.invert_config, invert_out, self.threads)

    def invert_labelmap(self):

        if not self.config.get('label_map_path'):
            logging.info('no label map found in config')
            return

        labelmap = join(self.proj_dir, self.config['label_map_path'])
        if not os.path.isfile(labelmap):
            logging.info('labelmap: {} not found')
            return

        label_inversion_dir = join(self.outdir, self.config['inverted_labels'])

        ilm = InvertLabelMap(self.invert_config, labelmap, label_inversion_dir, threads=self.threads)
        ilm.run()
        final_inverted_lm_dir = ilm.last_invert_output_dir

        organ_names_file_name = self.config.get('organ_names')
        if organ_names_file_name:
            organ_names = join(self.proj_dir, self.config['organ_names'])
            if not os.path.isfile(organ_names):
                logging.info('could not find organ names file: {}. Organ volume calculation not calculated'.format(organ_names))
                return
        else:
            logging.info('organ_names not specified in config file. Organ volumes not calculated')
            return

        voxel_size = float(self.config.get('voxel_size'))

        organ_vol_outfile = join(self.outdir, ORGAN_VOLS_OUT)

        calculate_volumes(final_inverted_lm_dir, organ_names, organ_vol_outfile, voxel_size)

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

        mesh_dir = join(self.proj_dir, self.config['isosurface_dir'])
        if not os.path.isdir(mesh_dir):
            logging.info('Mesh directory: {} not found'.format(mesh_dir))

        iso_out = self.paths.make('inverted_isosurfaces')

        logging.info('mesh inversion started')
        for mesh_path in common.GetFilePaths(mesh_dir):
            im = InvertMeshes(self.invert_config, mesh_path, iso_out)
            im.run()

    def save_metadata(self, metadata_filename):
        metata_path = join(self.outdir, metadata_filename)
        with open(metata_path, 'w') as fh:
            fh.write(yaml.dump(self.out_metadata, default_flow_style=False))

    def run_registration_schedule(self, config):
        """
        Run the registrations specified in the config file
        :param config: Config dictionary
        :return: 0 for success, or an error message
        """
        if self.no_qc:
            self.qc_dir = None
        else:
            self.qc_dir = self.paths.make('qc')

        do_pairwise = True if self.config.get('pairwise_registration') else False

        stage_params = self.generate_elx_parameters(config, do_pairwise)

        # Make dir to put averages in
        avg_dir = self.paths.make('averages')

        # Folder to put registered images in. Do not force mkdir if restarting registration
        make = False if self.restart_stage else 'f'
        root_reg_dir = self.paths.make('root_reg_dir', mkdir=make)

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

        if not self.no_qc and not self.restart_stage:
            input_histogram_dir = self.paths.make('input_image_histograms', parent=self.qc_dir)
            make_histograms(moving_vols_dir, input_histogram_dir)

        if not do_pairwise:
            fixed_vol = os.path.join(self.proj_dir, config.get('fixed_volume'))

        # Create a folder to store mid section coronal images to keep an eye on registration process

        if not self.no_qc:
            qc_image_dir = self.paths.make('qc_registered_images', parent=self.qc_dir)
            qc_metric_dir = self.paths.make('metric_charts', parent=self.qc_dir)

        reg_stages = config['registration_stage_params']
        reg_stage_ids = [s['stage_id'] for s in reg_stages]

        if self.restart_stage:
            if self.restart_stage not in reg_stage_ids:
                logging.error('It was requested to start at stage "{}", but this stage was not found in the config file')
                sys.exit()
            else:
                for i, s in enumerate(reg_stages):
                    if s['stage_id'] == self.restart_stage:
                        # Get the target and moving images from the previous stage
                        logging.info('restarting a previous registration from stage {}'.format(self.restart_stage))
                        fixed_vol, moving_vols_dir = self.get_volumes_for_restart(stage_params, self.restart_stage)
                        # Trim the previous stages
                        reg_stages = reg_stages[i:]

        if do_pairwise:
            logging.info('Using do_pairwise registration')
            RegMethod = PairwiseBasedRegistration
        else:
            logging.info('using target-based registration')
            RegMethod = TargetBasedRegistration

        for i, reg_stage in enumerate(reg_stages):

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
            if i == 0 and not self.restart_stage:
                fixed_mask = self.config.get('fixed_mask')

            # If we are doing target-based phenotype detection, we can used the fixed mask for every stage
            elif not self.config.get('generate_new_target_each_stage'):
                fixed_mask = self.config.get('fixed_mask')
            else:
                fixed_mask = None

            if fixed_mask:
                fixed_mask = join(self.config_dir, fixed_mask)

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
            if not do_pairwise:
                registrator.set_target(fixed_vol)
            registrator.run()
            # Make average
            average_path = join(avg_dir, '{0}.{1}'.format(stage_id, filetype))
            registrator.make_average(average_path)

            if not self.no_qc:
                stage_qc_image_dir = self.paths.make(join(qc_image_dir, stage_id))
                self.make_qc_images(stage_dir, stage_qc_image_dir)

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

        # Normalise linearly to the mean of the rois
        if config.get('normalisation_roi'):

            # Pass the final reg stage to be normalised
            norm_roi = config.get('normalisation_roi')
            #
            # self.normalise_registered_images(stage_dir, norm_dir, norm_roi)
            # if not self.no_qc:
            #     norm_qc_out_dir = self.paths.make('roi_region_overlays', parent=self.qc_dir)
            #     make_normalization_roi_qc_images(norm_dir, norm_roi, norm_qc_out_dir)
            #
            #     # Produce histograms of normalised images
            #     registered_normalised_hist_dir = self.paths.make('registered_normalised_histograms', parent=self.qc_dir)
            #     make_histograms(norm_dir, registered_normalised_hist_dir)

        if config.get('generate_deformation_fields'):
            make_vectors = not config.get('skip_deformation_fields')
            make_deformations_at_different_scales(config, root_reg_dir, self.outdir, make_vectors, self.threads,
                                                  filetype=config.get('filetype'), skip_histograms=config.get('no_qc'))

        if config.get('glcms'):
            self.create_glcms()

        logging.info("### Registration finished ###")

    def get_volumes_for_restart(self, stage_params, restart_stage):
        """
        If registration has restarted at a later stage, get the correct fixed and moving volumes

        Parameters
        ----------
        reg_stages: list
            A bunch of dicts, each specifying a registration stage
        restart_stage: str
            the stage_id the restart at
        """

        previous_id = None
        p = None
        for key, stage in stage_params.items():
            if key == restart_stage:
                previous_id = p
                break
            else:
                p = key

        if not previous_id:
            logging.error("'restart_at_stage:' refers to a stage that does not seem to exist")
            sys.exit()

        if self.config.get('generate_new_target_each_stage'):
            target = join(self.paths.get('averages'), previous_id + '.' + self.config['filetype'])
        else:
            if self.config.get('pad_dims'):
                target_dir = self.paths.get('padded_target')
                fixed_basename = splitext(basename(self.config.get('fixed_volume')))[0]
                target = join(target_dir, fixed_basename + '.' + self.config['filetype'])
            else:
                target = join(self.proj_dir, self.config.get('fixed_volume'))

        moving_vols_dir = join(self.paths.get('root_reg_dir'), previous_id)

        return target, moving_vols_dir

    def make_qc_images(self, im_dir, out_dir):
        """
        Generate a mid section slice to keep an eye on the registration process
        """
        for img_path in common.GetFilePaths(im_dir):
            loader = common.LoadImage(img_path)
            if not loader:
                logging.error('error making qc image: {}'.format(loader.error_msg))
            cast_img = sitk.Cast(sitk.RescaleIntensity(loader.img), sitk.sitkUInt8)
            arr = sitk.GetArrayFromImage(cast_img)
            slice_ = np.flipud(arr[:, :, arr.shape[2]/2])
            out_img = sitk.GetImageFromArray(slice_)
            base = splitext(basename(img_path))[0]
            out_path = join(out_dir, base + '.png')
            sitk.WriteImage(out_img, out_path, True)

    def create_glcms(self):
        """
        Create grey level co-occurence matrices. This is done in the main registration pipeline as we don't
        want to have to create GLCMs for the wildtypes multiple times when doing phenotype detection
        """
        if not glcm3d:
            return
        glcm_out_dir = self.paths.make('glcms')  # The vols to create glcms from
        registered_output_dir = join(self.outdir, self.config['normalised_output'])
        glcm3d.itk_glcm_generation(registered_output_dir, glcm_out_dir)

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
                        for k, v in reg_stage['elastix_parameters'].iteritems():
                            inherited_params[k] = v
                        # Now replace the stage-specific params with the merged new param dict
                        reg_stage['elastix_parameters'] = inherited_params

                if not referenced_stage_found:
                    sys.exit("Cannot inherit parameters from the stage: {}. Stage Not found. Check the config file".format(
                        stage_params['inherit_elx_params']))

            # Merge the global and stage-specific parameters
            global_params = config['global_elastix_params']
            if do_pairwise:
                global_params['WriteResultImage'] = 'false'
            param_vals.extend([global_params, reg_stage['elastix_parameters']])

            for p in param_vals:
                for param_name, param_value in p.iteritems():
                    if isinstance(param_value, list):  # parameter with multiple values for each resolution
                        # check whether we have didgets or strings, the latter need quoting
                        if all(are_numbers(v) for v in param_value):
                            val_list = [str(s).format() for s in param_value]
                            val_string = ' '.join(val_list)
                        else:
                            val_list = [str(s).format() for s in param_value]
                            val_quoted = ['"{}"'.format(s) for s in val_list]
                            val_string = ' '.join(val_quoted)
                        elxparams_formated.append('({0}  {1})\n'.format(param_name, val_string))
                    else:  # Find a nicer way to do this
                        if are_numbers(param_value):
                            elxparams_formated.append('({0}  {1})\n'.format(param_name, param_value))
                        else:
                            elxparams_formated.append('({0}  "{1}")\n'.format(param_name, param_value))

            stage_params[stage_id] = ''.join(elxparams_formated)
        return stage_params

    @staticmethod
    def test_elastix_installation():
        try:
            subprocess.check_output(['elastix'])
        except Exception:  # can't seem to log CalledProcessError
            logging.error('It looks like elastix may not be installed on your system')
            sys.exit()

    def get_config(self):
        return self.config

    def pad_inputs_and_modify_config(self):
        """
        Pad the input volumes, masks and labels
        if config['pad_dims'] == iterable: padd to these dimensions
        else: pad to the largest dimensions amongst the inputs.

        If padding occurs, update the config to point to the padded
        """
        replacements = {}

        config = self.config
        filetype = config.get('filetype')

        input_vol_dir = join(self.proj_dir, config['inputs'])

        fixed_vol_path = join(self.proj_dir, config['fixed_volume'])

        if os.path.isdir(input_vol_dir):
            input_vol_paths = common.GetFilePaths(input_vol_dir)
        else:
            input_vol_paths = common.get_inputs_from_file_list(input_vol_dir, self.config_dir)

        all_image_paths = [fixed_vol_path] + input_vol_paths

        try:
            iter(config['pad_dims'])  # Is it a list of dim lengths?
        except TypeError:
            maxdims = find_largest_dim_extents(all_image_paths)
        else:
            maxdims = config['pad_dims']

        replacements['pad_dims'] = str(maxdims)

        # pad the moving vols
        padded_mov_dir = self.paths.make('padded_inputs', 'f')

        # Pad the moving volumes and get a dict report of the padding occured
        moving_pad_info = pad_volumes(input_vol_paths, maxdims, padded_mov_dir, filetype)
        self.add_full_res_paths(moving_pad_info)

        config['inputs'] = padded_mov_dir
        replacements['inputs'] = relpath(padded_mov_dir, self.config_dir)

        # Create dir to put in padded volumes related to target such as mask and labelmap
        padded_fixed_dir = self.paths.make('padded_target', 'f')

        # Pad the fixed vol
        fixed_vol = self.config['fixed_volume']
        fixed_basename = splitext(basename(fixed_vol))[0]
        padded_fixed = join(padded_fixed_dir, '{}.{}'.format(fixed_basename, filetype))
        config['fixed_volume'] = padded_fixed
        fixed_vol_abs = join(self.proj_dir, fixed_vol)
        pad_volumes([fixed_vol_abs], maxdims, padded_fixed_dir, filetype)
        replacements['fixed_volume'] = relpath(padded_fixed, self.config_dir)

        # Pad labelmap, if pesent
        if config.get('label_map_path'):
            labels_path = config['label_map_path']
            label_basename = splitext(basename(labels_path))[0]
            padded_labels = join(padded_fixed_dir,'{}.{}'.format(label_basename, filetype))
            config['label_map_path'] = padded_labels
            labels_abs = join(self.proj_dir, labels_path )
            pad_volumes([labels_abs], maxdims, padded_fixed_dir, filetype)
            replacements['label_map_path'] = relpath(padded_labels, self.config_dir)
        else:
            logging.info("No labelmap path specified")

        if config.get('fixed_mask'):
            mask_path = config['fixed_mask']
            mask_basename = splitext(basename(mask_path))[0]
            padded_mask = join(padded_fixed_dir, '{}.{}'.format(mask_basename, filetype))
            config['fixed_mask'] = padded_mask
            mask_abs = join(self.proj_dir, mask_path)
            pad_volumes([mask_abs], maxdims, padded_fixed_dir, filetype)
            replacements['fixed_mask'] = relpath(padded_mask, self.config_dir)
        else:
            logging.info("No fixed mask specified")

        # If normalisation coordinates present, change coordinates appropriately
        # Need to find how many pixels have been added relative to target volume size

        norm_roi = config.get('normalisation_roi')
        if norm_roi:
            roi_starts = norm_roi[0]
            roi_ends = norm_roi[1]
            # TODO. this code is replicated im padd_volumes. Combine it
            target = sitk.ReadImage(fixed_vol_path)
            target_dims = target.GetSize()
            # The voxel differences between the vol dims and the max dims
            diffs = [m - v for m, v in zip(maxdims, target_dims)]

            # How many pixels to add to the upper bounds of each dimension, divide by two and round up to nearest int
            upper_extend = [d / 2 for d in diffs] # what about neagative values?

            # In case of differnces that cannot be /2. Get the remainder to add to the lower bound
            remainders = [d % 2 for d in diffs]

            # Add the remainders to the upper bound extension to get the lower bound extension
            le = [u + r for u, r in zip(upper_extend, remainders)]
            lower_extend = [l if l > 0 else 0 for l in le]
            new_roi_starts = [s + ex for s, ex in zip(roi_starts, lower_extend)]
            new_roi_ends = [e + ex for e, ex in zip(roi_ends, lower_extend)]
            config['normalisation_roi'] = [new_roi_starts, new_roi_ends]
            replacements['normalisation_roi'] = [new_roi_starts, new_roi_ends]

        if self.create_modified_config:
            config_name, ext = os.path.splitext(self.config_path)
            modified_config_path = "{}_pheno_detect{}".format(config_name, ext)
            replace_config_lines(self.config_path, replacements, modified_config_path)

    def add_full_res_paths(self, pad_info):
        """
        Generate a yaml file that contains the amount of padding applied to each image
        Also add the path to the full resolution image for each volume
        This file is used for inverting rois back onto the unpadded, full resolution images

        Parameters
        ----------
        pad_info: dict
            vol_id(minus extension):[(padding top xyz), (padding bottom xyz)]
        -------

        """

        full_res_root_folder = self.config.get('full_resolution_folder')
        if not full_res_root_folder:
            return

        full_res_subfolders = os.listdir(full_res_root_folder)
        if len(full_res_subfolders) < 1:
            return

        subfolder = self.config.get('full_resolution_subfolder_name')

        for vol_id in pad_info['data']:
            for f in full_res_subfolders:
                if f in vol_id:
                    if subfolder:
                        path = join(full_res_root_folder, f, subfolder)
                    else:
                        path = join(full_res_root_folder, f)
                    # Add path if it exists
                    if os.path.isdir(path):
                        pad_info['data'][vol_id]['full_res_folder'] = f
                    else:
                        pad_info['data'][vol_id]['full_res_folder'] = None
                        logging.warn("Cannot find full resolution path: {}".format(path))

        pad_info['root_folder'] = full_res_root_folder
        pad_info['full_res_subfolder'] = subfolder
        pad_info['log_file_endswith'] = self.config.get('log_file_endswith')
        if self.config.get('full_res_voxel_size'):
            pad_info['full_res_voxel_size'] = self.config.get('full_res_voxel_size')
        elif self.config.get('voxel_size_log_entry'):
            pad_info['voxel_size_log_entry'] = self.config.get('voxel_size_log_entry')
        pad_info_out = join(self.outdir, PAD_INFO_FILE)

        with open(pad_info_out, 'w') as fh:
            fh.write(yaml.dump(pad_info.to_dict()))

def replace_config_lines(config_path, key_values, config_path_out=None):
    """
    Replace lines in the config file. Did this rather than just writing out the config as we can't specify the
    order of elements from a dict and we lose comments too.

    This is also used by phenodetect.py
    """
    keys = key_values.keys()
    lines = []
    with open(config_path) as yif:
        for line in yif.readlines():
            for key in keys:
                if line.startswith(key):
                    line = '{}: {}\n'.format(key, key_values[key])
                    break
            lines.append(line)

    if config_path_out: # We rename the config
        c_out = config_path_out
    else:
        c_out = config_path

    with open(c_out, 'w') as yof:
        for outline in lines:
            yof.write(outline)


def set_origins_and_spacing(volpaths):
    for vol in volpaths:
        im = sitk.ReadImage(vol)
        if im.GetSpacing() != SPACING or im.GetOrigin() != ORIGIN:
            im.SetSpacing(SPACING)
            im.SetOrigin(ORIGIN)
            sitk.WriteImage(im, vol)


def parse_yaml_config(configfile):

    try:
        config = yaml.load(open(configfile, 'r'))
    except Exception as e:
        sys.exit("can't read the YAML config file - {}".format(e))
    return config


def are_numbers(value):
    """
    Is 'value' and int of a float. Used in formatting the elastix parameters
    :param value:
    :return:
    """

    try:
        float(value)
    except ValueError:
        pass
    else:
        return True

    try:
        int(value)
    except ValueError:
        pass
    else:
        return True

    return False


def find_largest_dim_extents(volpaths):
    merged = list(itertools.chain(volpaths))
    max_dims = None

    for path in merged:
        im = sitk.ReadImage(path)
        dims = im.GetSize()
        if not max_dims:
            max_dims = dims
        else:
            max_dims = [max(d[0], d[1]) for d in zip(dims, max_dims)]

    return max_dims


def make_histograms(in_dir, out_dir):
    """
    Make histograms of a series of input images and output to a html file
    """
    hist_batch(in_dir, out_dir, remove_zeros=True)


def mkdir_force(dir_):
    if os.path.isdir(dir_):
        shutil.rmtree(dir_)
    os.mkdir(dir_)


def mkdir_if_not_exists(dir_):
    if not os.path.exists(dir_):
        os.makedirs(dir_)

if __name__ == "__main__":


    parser = argparse.ArgumentParser("The MRC Harwell image registration pipeline")
    parser.add_argument('-c', dest='config', help='Config file (YAML format)', required=True)
    args = parser.parse_args()


    RegistraionPipeline(args.config)

