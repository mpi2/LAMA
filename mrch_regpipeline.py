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
    output_metadata_file: output_metadata.yaml  # this file is used by pheno_detect.py
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
from os.path import join, splitext, basename
import subprocess
import shutil
import sys
import argparse
import copy
import itertools
import logging
import datetime

import SimpleITK as sitk

import harwellimglib as hil
import yaml
from normalise import normalise
from invert import BatchInvertLabelMap, batch_invert_transform_parameters
import common


LOG_FILE = 'registration.log'
ELX_PARAM_PREFIX = 'elastix_params_'               # Prefix the generated elastix parameter files
INDV_REG_METADATA = 'reg_metadata.yaml'            # file name  for singleregistration metadata file
VOLUME_CALCULATIONS_PATH = 'organ_volumes.csv'
INVERT_ELX_TFORM_DIR = 'inverted_elx_tforms'

# Set the spacing and origins before registration
SPACING = (1.0, 1.0, 1.0)
ORIGIN = (0.0, 0.0, 0.0)


class RegistraionPipeline(object):
    def __init__(self, configfile, phenotyping=False):
        """This is the main function that is called by the GUI or from the command line.
        Reads in the config file, Creates directories, and initialises the registration process

        :param config: either a yaml file from the command line. Or a dict from the GUI
        :param: callback, function to return messages back to the GUI
        :return: nothing yet
        """
        self.proj_dir = os.path.dirname(configfile)
        config = parse_yaml_config(configfile)

        self.filetype = config['filetype']
        self.threads = str(config['threads'])
        self.voxel_size = float(config['voxel_size'])


        self.out_metadata = {'reg_stages': []}



        # all paths are relative to the config file directory
        self.config_dir = os.path.split(os.path.abspath(configfile))[0]
        self.outdir = join(self.config_dir, config['output_dir'])
        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)
        self.label_inversion_dir = join(self.outdir, 'label_inversion')
        self.add_metadata_path(self.label_inversion_dir, 'label_inversion_dir')

        self.volume_calculations_path = join(self.outdir, VOLUME_CALCULATIONS_PATH)
        logpath = join(self.outdir, LOG_FILE)
        common.init_log(logpath, 'Harwell registration pipeline')
        common.log_time("Registration started")

        inputvols_dir = join(self.config_dir, config['inputvolumes_dir'])

        if phenotyping:  # need to get fixed from wt project dir
            fixed_vol = join(config['wt_proj_dir'], config['fixed_volume'])
            if config.get('organ_names'):
                self.organ_names = join(config['wt_proj_dir'], config['organ_names'])
            else:
                self.organ_names = None
        else:
            fixed_vol = join(self.config_dir, config['fixed_volume'])
            if config.get('organ_names'):
                self.organ_names = join(self.config_dir, config['organ_names'])
            else:
                self.organ_names = None


        average_dir = join(self.outdir, 'averages')
        config['average_dir'] = average_dir
        mkdir_force(average_dir)

        # Pad out the moving dimensions of all the volumes to the same size. Do masks as well if required
        self.pad_dims = config.get('pad_dims')

        input_vol_paths = hil.GetFilePaths(inputvols_dir)
        input_vol_paths.append(fixed_vol)

        try:
            iter(self.pad_dims)
        except TypeError:
            maxdims = find_largest_dim_extents(input_vol_paths)
        else:
            maxdims = self.pad_dims

        # Modify the config file to add the padded paths or specified vols to the first stage.
        first_stage_config = config['registration_stage_params'][0]

        if self.pad_dims:
            # pad the moving vols
            paddedvols_dir = join(self.outdir, 'padded_inputs')
            mkdir_force(paddedvols_dir)
            pad_volumes(inputvols_dir, maxdims, paddedvols_dir, self.filetype)
            fixed_vol_dir = os.path.split(fixed_vol)[0]
            padded_fixed_dir = join(self.outdir, 'padded_target')
            mkdir_force(padded_fixed_dir)
            fixed_basename = splitext(basename(fixed_vol))[0]
            padded_fixed = join(padded_fixed_dir, '{}.{}'.format(fixed_basename, self.filetype))
            first_stage_config['fixed_vol'] = padded_fixed
            first_stage_config['movingvols_dir'] = paddedvols_dir
            pad_volumes(fixed_vol_dir, maxdims, padded_fixed_dir, self.filetype)
            self.add_metadata_path(padded_fixed, 'fixed_volume')

            if config.get('label_map'):
                lm = config['label_map'].get('path')
                if lm:
                    label_basename = splitext(basename(config['label_map'].get('path')))[0]
                    padded_labels = join(padded_fixed_dir,'{}.{}'.format(label_basename, self.filetype))
                    config['label_map']['path'] = padded_labels
            if config.get('fixed_mask'):
                mask_basename = splitext(basename(config['fixed_mask']))[0]
                padded_mask = join(padded_fixed_dir, '{}.{}'.format(mask_basename, self.filetype))
                config['fixed_mask'] = padded_mask
        else:
            first_stage_config['fixed_vol'] = fixed_vol
            first_stage_config['movingvols_dir'] = inputvols_dir
            self.add_metadata_path(fixed_vol, 'fixed_volume')
            config['label_map']['path'] = join(self.proj_dir,  config['label_map']['path'])

        self.add_metadata_path(config.get('fixed_mask'), 'fixed_mask')
        self.out_metadata['pad_dims'] = maxdims

        self.run_registration_schedule(config)

        metadata_filename = join(self.outdir, config['output_metadata_file'])

        tform_invert_dir = join(self.outdir, INVERT_ELX_TFORM_DIR)
        self.add_metadata_path(tform_invert_dir, 'inverted_elx_dir')
        self.save_metadata(metadata_filename)

        self._invert_elx_transform_parameters(metadata_filename, tform_invert_dir)

        if config.get('label_map'):
            self.invert_labelmap(config, tform_invert_dir)


    def _invert_elx_transform_parameters(self, metadata_filename, invert_out):
        """
        Invert the elastix output transform parameters. The inverted parameter files can then be used for inverting
        labelmaps and statistics overlays etc.
        """
        batch_invert_transform_parameters(metadata_filename, invert_out)

    def invert_labelmap(self, config, tform_invert_dir):

        labelmap = config['label_map'].get('path')

        if not labelmap:
            logging.warn('label_map path not present in config file')
            return

        invert_config = join(tform_invert_dir, 'invert.yaml')

        BatchInvertLabelMap(invert_config, labelmap, self.label_inversion_dir, self.organ_names,
                            do_organ_vol_calcs=True, threads=self.threads)


    def save_metadata(self, metadata_filename):
        metata_path = join(self.outdir, metadata_filename)
        with open(metata_path, 'w') as fh:
            fh.write(yaml.dump(self.out_metadata, default_flow_style=False))

    def add_metadata_path(self, path, key):
        """
        add paths of project directories and files to the output metadata dict
        :param path:
        :return:
        """
        if path:
            path = os.path.relpath(path, self.outdir)
        self.out_metadata[key] = path

    def validate_config(self, config):
        """
        Do some checks on the config file to check for errors
        :param config:
        :return: list, empty or containing errors
        """
        report = []
        required_params = ['output_dir', 'fixed_volume',
                           'inputvolumes_dir', 'filetype',
                           'global_elastix_params', 'registration_stage_params']

        for p in required_params:
            if p not in config:
                report.append("Entry: {} is required in the config file".format(p))

        stages = config['registration_stage_params']

        if len(stages) < 1:
            report.append("No stages specified")

        # ### Check for correct paths
        # if not os.path.isdir(config['inputvolumes_dir']):
        #     report.append("Input directory does not exit")
        #
        # if not os.path.isdir(config['output_dir']):
        #     report.append("Output directory does not exit")

        return report

    def run_registration_schedule(self, config):
        """
        Run the registrations specified in the config file
        :param config: Config dictionary
        :return: 0 for success, or an error message
        """

        # Testing: reapeat the registration stage replacing the target with the newly generated average.
        repeat_registration = True

        regenerate_target = config['generate_new_target_each_stage']

        for i, reg_stage in enumerate(config['registration_stage_params']):

            #  Make the stage output dir
            stage_id = reg_stage['stage_id']
            self.out_metadata['reg_stages'].append(stage_id)
            stage_dir = join(self.outdir, stage_id)
            mkdir_force(stage_dir)

            print("### Current registration step: {} ###".format(stage_id))

            # Make the elastix parameter file for this stage
            elxparam = self.make_elastix_parameter_file(config, reg_stage)
            elxparam_path = join(stage_dir, ELX_PARAM_PREFIX + stage_id + '.txt')

            with open(elxparam_path, 'w') as fh:
                if elxparam:  # Not sure why I put this here
                    fh.write(elxparam)

            # For rigid, miss out registering to self

            if stage_id == 'rigid':
                fixed_mask = config.get('fixed_mask')

            # Do the registration
            self.elx_registration(elxparam_path,
                                  reg_stage['fixed_vol'],
                                  reg_stage['movingvols_dir'],
                                  stage_dir,
                                  fixed_mask
                                  )

            # Normalise, if required
            if reg_stage.get('normalise_registered_output'):
                self.normalise_registered_images(reg_stage, stage_id, stage_dir)

            # Generate deformation fields
            if reg_stage.get('do_analysis'):
               self.do_analysis(stage_dir)

            # Make average
            average_path = join(config['average_dir'], '{0}.{1}'.format(stage_id, self.filetype))
            make_average(stage_dir, average_path)

            # Create new avergae
            if config.get('re-register_each_stage'):
                average_path, stage_dir = self.repeat_registration(average_path, reg_stage['movingvols_dir'], elxparam_path, average_path)

            # Setup the fixed, moving and mask_dir for the next stage, if there is one
            if i + 1 < len(config['registration_stage_params']):
                next_reg_stage = config['registration_stage_params'][i + 1]

                if regenerate_target:
                    next_reg_stage['fixed_vol'] = average_path  # The avergae from the previous step
                else:
                    next_reg_stage['fixed_vol'] = reg_stage['fixed_vol']  # The same fixed volume as thje previous step

                next_reg_stage['movingvols_dir'] = stage_dir  # The output dir of the previous registration

        print("### Registration finished ###")
        common.log_time('gistration Finished')

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

        return (new_average_path, output_dir)


    def do_analysis(self, stage_dir):

        deformation_dir = join(self.outdir, 'deformation_fields')
        jacobians_dir = join(self.outdir, 'spatial_jacobians')
        mkdir_force(deformation_dir)
        mkdir_force(jacobians_dir)
        self.generate_deformation_fields(stage_dir, deformation_dir, jacobians_dir)
        self.add_metadata_path(deformation_dir, 'deformation_fields')
        self.add_metadata_path(jacobians_dir, 'jacobians')

    def normalise_registered_images(self, reg_stage_config, stage_id, stage_dir):

        norm_settings = reg_stage_config.get('normalise_registered_output')
        roi_starts = norm_settings[0]
        roi_ends = norm_settings[1]
        norm_dir = join(self.outdir, stage_id + "_" + 'normalised')
        mkdir_force(norm_dir)

        logging.info('Normalised registered output using ROI: {}:{}'.format(
            ','.join([str(x) for x in roi_starts]), ','.join([str(x) for x in roi_ends])))
        normalise(stage_dir, norm_dir, roi_starts, roi_ends)

        self.add_metadata_path(norm_dir, 'normalized_registered')

    def generate_deformation_fields(self, registration_dir, deformation_dir, jacobian_dir):
        """
        Run transformix on the specified registration stage to generate deformation fields and spatial jacobians
        :param registration_dir:
        :param deformation_dir:
        :param jacobian_dir:
        :return:
        """
        print('### Generating deformation files ###')
        print registration_dir, deformation_dir, jacobian_dir
        # Get the transform parameters
        dirs_ = os.listdir(registration_dir)
        print dirs_
        for dir_ in dirs_:
            if os.path.isdir(join(registration_dir, dir_)):

                elx_tform_file = join(registration_dir, dir_, 'TransformParameters.0.txt')
                if not os.path.isfile(elx_tform_file):
                    sys.exit("### Error. Cannot find elastix transform parameter file: {}".format(elx_tform_file))

                cmd = ['transformix',
                       '-tp', elx_tform_file,
                       '-out', deformation_dir,
                       '-def', 'all',
                       '-jac', 'all']

                try:
                    output = subprocess.check_output(cmd)
                except Exception as e:  # Can't seem to log CalledProcessError
                    logging.warn('transformix failed {}'.format(', '.join(cmd)))
                    sys.exit('### Transformix failed ###\nError message: {}\nelastix command:{}'.format(e, cmd))
                else:
                    # Rename the outputs and move to correct dirs
                    volid = basename(dir_)
                    deformation_out = join(deformation_dir, 'deformationField.{}'.format(self.filetype))
                    jacobian_out = join(deformation_dir, 'spatialJacobian.{}'.format(self.filetype))
                    deformation_renamed = join(deformation_dir,
                                                       '{}_deformationFied.{}'.format(volid, self.filetype))
                    jacobian_renamed = join(jacobian_dir, '{}_spatialJacobian.{}'.format(volid, self.filetype))
                    shutil.move(deformation_out, deformation_renamed)
                    shutil.move(jacobian_out, jacobian_renamed)

    def make_elastix_parameter_file(self, config, stage_params):
        """
        :param main_elx_params: Dict with elastix params for all stages
        :param stage_params: Dict with stage-specific elastix params
        :return: str, formatted elastix parameters
        """
        elxparams_formated = []
        param_vals = []  # The parameters to format

        if 'inherit_elx_params' in stage_params:  # Use another stage's params as a starting point

            # Find the referenced stage's parameters
            referenced_stage_found = False
            for stage in config['registration_stage_params']:
                if stage['stage_id'] == stage_params['inherit_elx_params']:
                    referenced_stage_found = True
                    inherited_params = copy.deepcopy(stage['elastix_parameters'])
                    #  Overwrite the inherited params with this stage-specific ones
                    for k, v in stage_params['elastix_parameters'].iteritems():
                        inherited_params[k] = v
                    # Now replace the stage-specific params with the merged new param dict
                    stage_params['elastix_parameters'] = inherited_params

            if not referenced_stage_found:
                sys.exit("Cannot inherit parameters from the stage: {}. Stage Not found. Check the config file".format(
                    stage_params['inherit_elx_params']))

        # Merge the global and stage-specific parameters
        param_vals.extend([config['global_elastix_params'], stage_params['elastix_parameters']])

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

        return ''.join(elxparams_formated)

    def elx_registration(self, elxparam_file, fixed, movdir, stagedir, fixed_mask=None):
        """
        Run elastix for one stage of the registration pipeline

        """
        movlist = hil.GetFilePaths(movdir)

        for mov in movlist:
            mov_basename = splitext(basename(mov))[0]
            outdir = join(stagedir, mov_basename)
            mkdir_force(outdir)

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
            except Exception:  # can't seem to log CalledProcessError
                logging.exception('registration falied:\n\n{}\n\n'.format(cmd))
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

                # Apply the transform to the mask so it can be used in the next stage
                tp_file = join(outdir, "TransformParameters.0.txt")
                # transform_mask(moving_mask, mask_dir_out, tp_file)

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



def pad_volumes(voldir, max_dims, outdir, filetype):
    """
    :param maxdims: list, dimensions to pad to (X, Y, Z)
    All volumes should be the same dimensions to aid in downstraem analysis and to prevent clipping when moving is
    larger than the fixed
    :return:
    """

    volpaths = hil.GetFilePaths(voldir)

    if len(volpaths) < 1:
        sys.exit("Can't find any volumes in {}".format(voldir))

    print 'Padding to {} - {} volumes/masks:'.format(str(max_dims), str(len(volpaths)))

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
    avg = hil.Average(moving_mask_dir)
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
    vols = hil.GetFilePaths(vol_dir)

    average = hil.Average(vols)

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser("The MRC Harwell image registration pipeline")
    parser.add_argument('-c', dest='config', help='Config file (YAML format)', required=True)
    args = parser.parse_args()

    RegistraionPipeline(args.config)
