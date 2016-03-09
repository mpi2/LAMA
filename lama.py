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
    inputvolumes_dir: inputs  # directory with moving volumes
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
from collections import OrderedDict
import copy
import itertools
import logging
import yaml
import SimpleITK as sitk
import numpy as np
from normalise import normalise
from invert import InvertLabelMap, InvertMeshes, batch_invert_transform_parameters
import common
try:
    import glcm3d
except ImportError:
    glcm3d = False
from calculate_organ_size import calculate_volumes
from validate_config import validate_reg_config
from deformations import generate_deformation_fields
from paths import RegPaths
from metric_charts import make_charts

LOG_FILE = 'LAMA.log'
ELX_PARAM_PREFIX = 'elastix_params_'               # Prefix the generated elastix parameter files
INDV_REG_METADATA = 'reg_metadata.yaml'            # file name  for singleregistration metadata file
INVERT_CONFIG = 'invert.yaml'
ORGAN_VOLS_OUT = 'organ_volumes.csv'



# Set the spacing and origins before registration
SPACING = (1.0, 1.0, 1.0)
ORIGIN = (0.0, 0.0, 0.0)


class RegistraionPipeline(object):
    def __init__(self, configfile, create_modified_config=True):
        """This is the main function that is called by the GUI or from the command line.
        Reads in the config file, Creates directories, and initialises the registration process

        :param config: either a yaml file from the command line. Or a dict from the GUI
        :param: callback, function to return messages back to the GUI
        :return: nothing yet
        """
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

        # Where to put all the output
        self.outdir = self.paths.make('outdir')

        # Number oif threads to use during elastix registration
        threads = self.config.get('threads')
        if threads:
            self.threads = str(threads)

        # The filtype extension to use for registration output
        self.filetype = config['filetype']

        logging.info(common.git_log())

        logging.info("Registration started")

        self.restart_stage = config.get('restart_at_stage')

        # Pad the inputs. Also changes the config object to point to these newly padded volumes

        if config.get('pad_dims') and not self.restart_stage:
            # so just always pad for now# : Fix
            self.pad_inputs_and_modify_config()

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

        batch_invert_transform_parameters(self.config_path, self.invert_config, invert_out)

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
            logging.info('Mesh directory: {} not found')

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

        stage_params = self.generate_elx_parameters(config)

        # Make dir to put averages in
        avg_dir = self.paths.make('averages')

        # Folder to put registered images in. Do not force mkdir if restarting registration
        make = False if self.restart_stage else 'f'
        root_reg_dir = self.paths.make('root_reg_dir', mkdir=make)

        # if True: create a new fixed volume by averaging outputs
        # if False: use the same fixed volume at each registration stage
        regenerate_target = config['generate_new_target_each_stage']
        if regenerate_target:
            logging.info('Creating new target each stage for population average creation')
        else:
            logging.info('Using same target for each stage')

        filetype = config['filetype']

        # Set the moving vol dir and the fixed image for the first satge
        moving_vols_dir = config['inputvolumes_dir']
        fixed_vol = os.path.join(self.proj_dir, config['fixed_volume'])

        # Create a folder to store mid section coroal images to keep an eye on registration process
        qc_dir = self.paths.make('qc')
        qc_image_dir = self.paths.make('qc_images', parent=qc_dir)
        qc_metric_dir = self.paths.make('metric_charts', parent=qc_dir)

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

            # For the first stage, we can use the fixedf mask for registration.
            # Sometimes helps with the 'too many samples map outside fixed image' problem
            if i == 0 and not self.restart_stage:
                fixed_mask = self.paths.get('fixed_mask')
            else:
                fixed_mask = None

            # Do the registrations
            self.elx_registration(elxparam_path,
                                  fixed_vol,
                                  moving_vols_dir,
                                  stage_dir,
                                  fixed_mask
                                  )

            stage_qc_image_dir = self.paths.make(join(qc_image_dir, stage_id))

            self.make_qc_images(stage_dir, stage_qc_image_dir)

            stage_metrics_dir = join(qc_metric_dir, stage_id)
            self.paths.make(stage_metrics_dir)
            make_charts(stage_dir, stage_metrics_dir)

            # Make average
            average_path = join(avg_dir, '{0}.{1}'.format(stage_id, filetype))
            make_average(stage_dir, average_path)

            # Reregister the inputs to this stage to the average from this stage to create a new average
            if config.get('re-register_each_stage'):
                average_path, stage_dir = self.repeat_registration(average_path, moving_vols_dir, elxparam_path, average_path)

            # Setup the fixed, moving and mask_dir for the next stage, if there is one
            if i + 1 < len(config['registration_stage_params']):
                if regenerate_target:
                    fixed_vol = average_path  # The avergae from the previous step

                moving_vols_dir = stage_dir  # The output dir of the previous registration

        # Normalise final output, if required
        if config.get('background_roi_zyx_norm'):
            self.normalise_registered_images(stage_dir, config.get('background_roi_zyx_norm')) # Pass the final reg stage to be normalised

        if config.get('generate_deformation_fields'):
            stage_to_start = config['generate_deformation_fields']
            def_stage_ids = reg_stage_ids[reg_stage_ids.index(stage_to_start):]
            def_stage_dirs = [join(root_reg_dir, x) for x in def_stage_ids]
            deformation_dir = self.paths.make('deformations')
            jacobians_dir = self.paths.make('jacobians')

            generate_deformation_fields(def_stage_dirs, deformation_dir, jacobians_dir)

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
            target = join(self.paths.get('averages'), previous_id)
        else:
            if self.config.get('pad_dims'):
                target_dir = self.paths.get('padded_target')
                fixed_basename = splitext(basename(self.config.get('fixed_volume')))[0]
                target = join(target_dir, fixed_basename + '.nrrd')
            else:
                target = join(self.proj_dir, self.config.get('fixed_volume'))

        moving_vols_dir = join(self.paths.get('root_reg_dir'), previous_id)

        return target, moving_vols_dir

    def make_qc_images(self, im_dir, out_dir):
        """
        Generate a mid section slice to keep an eye on the registration process
        """
        for img_path in common.GetFilePaths(im_dir):
            img = sitk.ReadImage(img_path)
            cast_img = sitk.Cast(sitk.RescaleIntensity(img), sitk.sitkUInt8)
            arr = sitk.GetArrayFromImage(cast_img)
            slice_ = np.flipud(arr[:, :, arr.shape[2]/2])
            out_img = sitk.GetImageFromArray(slice_)
            base = splitext(basename(img_path))[0]
            out_path = join(out_dir, base + '.png')
            sitk.WriteImage(out_img, out_path)

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

    def repeat_registration(self, population_average, input_dir, elxparam_path, average_path):
        """
        This is a temporary method for testing. For each registration stage, use the output average as the target.
        But instead of using the outputs of the stage as inputs, use the inputs from the stage. Only then use the
        outputs as targets and tp create an average for the next stage.

        Parameters
        ----------
        population_average: str
            path to target
        input_dir: str
            dir containg moving images

        Returns
        -------
        avergage: str
            path to tem
        output_dir: str
            path to registered images
        """
        output_dir = os.path.join(input_dir, 're-registered')
        mkdir_force(output_dir)

         # Do the registration
        self.elx_registration(elxparam_path,
                              population_average,
                              input_dir,
                              output_dir
                              )
        # Make average
        avg_path, ext = os.path.splitext(average_path)
        new_path = avg_path + 'REDONE'
        new_average_path = new_path + ext
        make_average(output_dir, new_average_path)

        return new_average_path, output_dir

    def normalise_registered_images(self, stage_dir, norm_roi):

        roi_starts = norm_roi[0]
        roi_ends = norm_roi[1]
        norm_dir = self.paths.make('intensity', mkdir='force')

        logging.info('Normalised registered output using ROI: {}:{}'.format(
            ','.join([str(x) for x in roi_starts]), ','.join([str(x) for x in roi_ends])))
        normalise(stage_dir, norm_dir, roi_starts, roi_ends)

    def generate_elx_parameters(self, config):
        """
        :param main_elx_params: Dict with elastix params for all stages
        :param stage_params: Dict with stage-specific elastix params
        :return: str, formatted elastix parameters
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
            param_vals.extend([config['global_elastix_params'], reg_stage['elastix_parameters']])

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

    def elx_registration(self, elxparam_file, fixed, movdir, stagedir, fixed_mask=None):
        """
        Run elastix for one stage of the registration pipeline

        """

        # If inputs_vols is a file get the specified root and paths from it
        if os.path.isdir(movdir):
            movlist = common.GetFilePaths(movdir)
        else:
            movlist = common.get_inputs_from_file_list(movdir)

        for mov in movlist:
            mov_basename = splitext(basename(mov))[0]
            outdir = self.paths.make(join(stagedir, mov_basename), 'f')

            set_origins_and_spacing([fixed, mov])

            cmd = ['elastix',
                   '-f', fixed,
                   '-m', mov,
                   '-out', outdir,
                   '-p', elxparam_file,
                   '-threads', self.threads]

            if fixed_mask:
                pass
                #cmd.extend(['-fMask', fixed_mask])

            try:
                subprocess.check_output(cmd)
            except Exception as e:  # can't seem to log CalledProcessError
                logging.exception('registration falied:\n\ncommand: {}\n\n error:{}'.format(cmd, e))
                raise
            else:
                # Rename the registered output.
                elx_outfile = join(outdir, 'result.0.{}'.format(self.filetype))
                new_out_name = join(outdir, '{}.{}'.format(mov_basename, self.filetype))
                shutil.move(elx_outfile, new_out_name)

                # add registration metadata
                reg_metadata_path = join(outdir, INDV_REG_METADATA)
                fixed_vol_relative = os.path.relpath(fixed, outdir)
                reg_metadata = {'fixed_vol': fixed_vol_relative}
                with open(reg_metadata_path, 'w') as fh:
                    fh.write(yaml.dump(reg_metadata, default_flow_style=False))

    def transform_mask(self, moving_mask, mask_dir_out, tp_file):
        """
        After a moving mask is used in registering a volume, it will be needed for the next stage. It needs to have the
         same transformation applied to it as the volume it is masking in order to still fit over the image correctly.
         Also compress the result as masks are highly compressible.
         Set interpolation mode to linear for binaray masks - (FinalBSplineInterpolationOrder 0)
        :return:
        """
        temp_tpfile = "_temp_mask_transform_file.txt"
        with open(tp_file, 'r') as tp_fh, open(temp_tpfile, 'w') as temp_fh:
            for line in tp_fh:
                if line.startswith('(FinalBSplineInterpolationOrder'):
                    line = '(FinalBSplineInterpolationOrder 0)\n'
                temp_fh.write(line)

        cmd = ['transformix',
               '-in', moving_mask,
               '-tp', temp_tpfile,
               '-out', mask_dir_out,
               ]

        try:
            subprocess.check_output(cmd)
        except Exception as e:  # Can't seem to log CalledProcessError
            logging.warn('transformix failed {}'.format(', '.join(cmd)))
            sys.exit('### Transformix failed while transforming mask ###\nelastix command:{}'.format(cmd))
        else:
            os.remove(temp_tpfile)
            # Rename the transformed mask
            outname = join(mask_dir_out, 'result.nrrd')
            newname = join(mask_dir_out, basename(moving_mask))
            mask = sitk.ReadImage(outname)
            sitk.WriteImage(mask, newname, True)
            os.remove(outname)

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

        input_vol_dir = join(self.proj_dir, config['inputvolumes_dir'])

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
        pad_volumes(input_vol_paths, maxdims, padded_mov_dir, filetype)
        config['inputvolumes_dir'] = padded_mov_dir
        replacements['inputvolumes_dir'] = relpath(padded_mov_dir, self.config_dir)

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

        norm_roi = config.get('background_roi_zyx_norm')
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
            config['background_roi_zyx_norm'] = [new_roi_starts, new_roi_ends]
            replacements['background_roi_zyx_norm'] = [new_roi_starts, new_roi_ends]

        if self.create_modified_config:
            config_name, ext = os.path.splitext(self.config_path)
            modified_config_path = "{}_pheno_detect{}".format(config_name, ext)
            replace_config_lines(self.config_path, replacements, modified_config_path)

def pad_roi(roi_starts, roi_ends):
    pass

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


def pad_volumes(volpaths, max_dims, outdir, filetype):
    """
    Pad volumes, masks, labels. Output files will have same name as original, but be in a new output folder

    Parameters
    ----------
    volpaths: list
        list of paths to volumes
    maxdims: list
        dimensions to pad to (z, Y, x)
    outdir: str
        path to output dir
    """
    logging.info('Padding to {} - {} volumes/masks:'.format(str(max_dims), str(len(volpaths))))

    for path in volpaths:

        vol = sitk.ReadImage(path)
        vol_dims = vol.GetSize()

        # The voxel differences between the vol dims and the max dims
        diffs = [m - v for m, v in zip(max_dims, vol_dims)]

        # How many pixels to add to the upper bounds of each dimension, divide by two and round up to nearest int
        upper_extend = [d / 2 for d in diffs]

        # In case of differnces that cannot be /2. Get the remainder to add to the lower bound
        remainders = [d % 2 for d in diffs]

        # Add the remainders to the upper bound extension to get the lower bound extension
        lower_extend = [u + r for u, r in zip(upper_extend, remainders)]

        # if any values are negative, stop. We need all volumes to be the same size
        for ex_val in zip(lower_extend, upper_extend):
            if ex_val[0] < 0 or ex_val[1] < 0:
                sys.exit("can't pad images. {} is larger than the specified volume size"
                         "Check the 'pad_dims' in the config file".format(path))

        # Pad the volume. New pixels set to zero
        padded_vol = sitk.ConstantPad(vol, upper_extend, lower_extend, 0)
        padded_vol.SetOrigin(ORIGIN)
        padded_vol.SetSpacing(SPACING)

        input_basename = splitext(basename(path))[0]
        padded_outname = join(outdir, '{}.{}'.format(input_basename, filetype))
        sitk.WriteImage(padded_vol, padded_outname, True)


def make_average_mask(moving_mask_dir, fixed_mask_path):
    """
    Create a fixed mask by creating a mask from moving masks. Just create an average, and binarize by setting all
    pixels that are not 0 to 1
    :param moving_mask_dir:
    :param fixed_mask_dir:
    :return:
    """
    avg = common.Average(moving_mask_dir)
    array = sitk.GetArrayFromImage(avg)
    array[array != 0] = 1
    avg_out = sitk.GetImageFromArray(array)
    sitk.WriteImage(avg_out, fixed_mask_path, True)


def set_origins_and_spacing(volpaths):
    return
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


def make_average(vol_dir, out_path):
    """
    Create an average of the the input embryo volumes.
    This will search subfolders for all the registered volumes within them
    """
    vols = common.GetFilePaths(vol_dir)

    average = common.Average(vols)

    sitk.WriteImage(average, out_path, True)  # Compressed=True
    # Check that it's been created
    if os.path.exists(out_path):
        return 0
    else:
        return 'Cannot make average at {}'.format(out_path)


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
