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



"""

import os
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

LOG_MODE = logging.DEBUG
LOG_FILE = 'registration.log'

class RegistraionPipeline(object):

    def __init__(self, configfile, phenotyping=False):
        """This is the main function that is called by the GUI or from the command line.
        Reads in the config file, Creates directories, and initialises the registration process

        :param config: either a yaml file from the command line. Or a dict from the GUI
        :param: callback, function to return messages back to the GUI
        :return: nothing yet
        """
        config = parse_yaml_config(configfile)

        config_failed = self.validate_config(config)
        if config_failed:
            sys.exit('\n\n'.join(config_failed))

        self.filetype = config['filetype']

        self.out_metadata = {}

        # all paths are relative to the config file directory
        self.config_dir = os.path.split(os.path.abspath(configfile))[0]
        self.outdir = os.path.join(self.config_dir, config['output_dir'])
        mkdir_force(self.outdir)

        logfile = os.path.join(self.outdir, LOG_FILE)
        logging.basicConfig(filename=logfile,level=LOG_MODE, filemode="w")
        logging.info("Harwell registration pipeline")
        self.log_time('Started')

        inputvols_dir = os.path.join(self.config_dir,  config['inputvolumes_dir'])

        if phenotyping:  # need to get fixed from wt project dir
            fixed_vol = os.path.join(config['wt_proj_dir'], config['fixed_volume'])

        else:
            fixed_vol = os.path.join(self.config_dir, config['fixed_volume'])

        average_dir = os.path.join(self.outdir, 'averages')
        config['average_dir'] = average_dir
        mkdir_force(average_dir)

        # Pad out the moving dimensions of all the volumes to the same size. Do masks as well if required
        pad_dims = config.get('pad_dims')

        # TODO: save either the given, or determined pad dims into output_metadata
        # todo: for pheno detection. always use the same dims as wt reg
        input_vol_paths = hil.GetFilePaths(inputvols_dir)
        input_vol_paths.append(fixed_vol)

        try:
            iter(pad_dims)
        except TypeError:
            maxdims = find_largest_dim_extents(input_vol_paths)
        else:
            maxdims = pad_dims

        if pad_dims:
            # pad the moving vols
            paddedvols_dir = os.path.join(self.outdir, 'padded_inputs')
            mkdir_force(paddedvols_dir)
            pad_volumes(inputvols_dir, maxdims, paddedvols_dir, self.filetype)

        # Modify the config file to add the padded paths or specified vols to the first stage.
        first_stage_config = config['registration_stage_params'][0]

        if pad_dims:
            fixed_vol_dir = os.path.split(fixed_vol)[0]
            padded_fixed_dir = os.path.join(self.outdir, 'padded_target')
            mkdir_force(padded_fixed_dir)
            basename = os.path.splitext(os.path.basename(fixed_vol))[0]
            mask_basename = os.path.splitext(os.path.basename(config['fixed_mask']))[0]
            padded_fixed = os.path.join(padded_fixed_dir, '{}.{}'.format(basename, self.filetype))
            padded_mask = os.path.join(padded_fixed_dir, '{}.{}'.format(mask_basename, self.filetype))
            config['fixed_mask'] = padded_mask
            first_stage_config['fixed_vol'] = padded_fixed
            first_stage_config['movingvols_dir'] = paddedvols_dir
            pad_volumes(fixed_vol_dir, maxdims, padded_fixed_dir, self.filetype)
            self.add_metadata_path(padded_fixed, 'fixed_volume')
        else:
            first_stage_config['fixed_vol'] = fixed_vol
            first_stage_config['movingvols_dir'] = inputvols_dir
            self.add_metadata_path(fixed_vol, 'fixed_volume')

        self.add_metadata_path(config['fixed_mask'], 'fixed_mask')
        self.out_metadata['pad_dims'] = maxdims

        self.do_registration(config)
        self.save_metadata(config['output_metadata_file'])

    def log_time(self, msg):
        now = datetime.datetime.now()
        logging.info("{}: {}/{}/{} - {}:{}".format(msg, now.day, now.month, now.year, now.hour, now.minute))

    def save_metadata(self, metadata_filename):
        metata_path = os.path.join(self.outdir, metadata_filename)
        with open(metata_path, 'w') as fh:
            fh.write(yaml.dump(self.out_metadata, default_flow_style=False))

    def add_metadata_path(self, path, key):
        """

        :param path:
        :return:
        """
        relpath = os.path.relpath(path, self.config_dir)
        self.out_metadata[key] = relpath

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

    def do_registration(self, config):
        """
        Run the registrations specified in the config file
        :param config: Config dictionary
        :return: 0 for success, or an error message
        """

        deformation_dir = os.path.join(self.outdir, 'deformation_fields')

        jacobians_dir = os.path.join(self.outdir, 'spatial_jacobians')
        regenerate_target = config['generate_new_target_each_stage']

        delete_stages = []
        for i, reg_stage in enumerate(config['registration_stage_params']):

            #  Make the stage output dir
            stage_id = reg_stage['stage_id']
            stage_dir = os.path.join(self.outdir, stage_id)
            mkdir_force(stage_dir)

            if reg_stage.get('delete_stage_files'):
                delete_stages.append((i + 1, stage_dir))

            print("### Current registration step: {} ###".format(stage_id))

            # Make the elastix parameter file for this stage
            elxparam = self.make_elastix_parameter_file(config, reg_stage)
            elxparam_path = os.path.join(stage_dir, 'elastix_params_' + stage_id + '.txt')

            with open(elxparam_path, 'w') as fh:
                if elxparam:  # Not sure why I put this here
                    fh.write(elxparam)

            # For rigid, miss out registering to self
            # Do a single stage registration

            if stage_id == 'rigid':
                fixed_mask = config.get('fixed_mask')

            self.elx_registration(elxparam_path,
                             reg_stage['fixed_vol'],
                             reg_stage['movingvols_dir'],
                             stage_dir,
                             config['threads'],
                             fixed_mask
                             )

            if reg_stage.get('normalise_registered_output'):

                norm_settings = reg_stage.get('normalise_registered_output')
                roi_starts = norm_settings[0]
                roi_ends = norm_settings[1]
                norm_dir = os.path.join(self.outdir, stage_id + "_" + 'normalised')
                mkdir_force(norm_dir)

                logging.info('Normalised registered output using ROI: {}:{}'.format(
                    ','.join([str(x) for x in roi_starts]), ': '.join([str(x) for x in roi_ends])))
                normalise(stage_dir, norm_dir, roi_starts, roi_ends)
                
                self.add_metadata_path(norm_dir, 'normalized_registered')

            # Generate deformation fields
            if reg_stage.get('do_analysis'):
                mkdir_force(deformation_dir)
                mkdir_force(jacobians_dir)
                self.generate_deformation_fields(stage_dir, deformation_dir, jacobians_dir)
                self.add_metadata_path(deformation_dir, 'deformation_fields')
                self.add_metadata_path(jacobians_dir, 'jacobians')

            # Make average
            average_path = os.path.join(config['average_dir'], '{0}.{1}'.format(stage_id, self.filetype))
            make_average(stage_dir, average_path)

            for x in delete_stages:
                if x[0] == i:
                    if os.path.isdir(x[1]):
                        shutil.rmtree(x[1])

            # Setup the fixed, moving and mask_dir for the next stage, if there is one
            if i + 1 < len(config['registration_stage_params']):
                next_reg_stage = config['registration_stage_params'][i + 1]

                if regenerate_target:
                    next_reg_stage['fixed_vol'] = average_path  # The avergae from the previous step
                else:
                    next_reg_stage['fixed_vol'] = reg_stage['fixed_vol']  # The same fixed volume as thje previous step

                if 'movingvols_dir' in next_reg_stage:  # Moving volumes are specified, don't just use previous stage
                    stage_output_to_use = next_reg_stage['movingvols_dir']
                    for stage in config['registration_stage_params']:
                        if stage['stage_id'] == stage_output_to_use:
                            next_reg_stage['movingvols_dir'] = os.path.join(self.outdir, stage_id)

                else:
                    next_reg_stage['movingvols_dir'] = stage_dir  # The output dir of the previous registration

        print("### Registration pipeline finished ###")
        self.log_time('Finished')

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
            if os.path.isdir(os.path.join(registration_dir, dir_)):

                elx_tform_file = os.path.join(registration_dir, dir_, 'TransformParameters.0.txt')
                if not os.path.isfile(elx_tform_file):
                    sys.exit("### Error. Cannot find elastix transform parameter file: {}".format(elx_tform_file))

                cmd = ['transformix',
                       '-tp', elx_tform_file,
                       '-out', deformation_dir,
                       '-def', 'all',
                       '-jac', 'all']

                try:
                    output = subprocess.check_output(cmd)
                except subprocess.CalledProcessError as e:
                    sys.exit('### Transformix failed ###\nError message: {}\nelastix command:{}'.format(output, cmd))
                else:
                    # Rename the outputs and move to correct dirs
                    volid = os.path.basename(dir_)
                    deformation_out = os.path.join(deformation_dir, 'deformationField.{}'.format(self.filetype))
                    jacobian_out = os.path.join(deformation_dir, 'spatialJacobian.{}'.format(self.filetype))
                    deformation_renamed = os.path.join(deformation_dir, '{}_deformationFied.{}'.format(volid, self.filetype))
                    jacobian_renamed = os.path.join(jacobian_dir, '{}_spatialJacobian.{}'.format(volid, self.filetype))
                    shutil.move(deformation_out, deformation_renamed)
                    shutil.move(jacobian_out, jacobian_renamed)

    def make_elastix_parameter_file(self, config, stage_params):
        """
        :param main_elx_params: Dict with elastix params for all stages
        :param stage_params: Dict with stage-specific elastix params
        :return: str, formatted elastix parameters
        """
        elxparams_formated = []
        param_vals = [] # The parameters to format

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

    def elx_registration(self, elxparam_file, fixed, movdir, stagedir, threads, fixed_mask=None):
        """
        Run elastix for one volume

        :return:
        """
        movlist = hil.GetFilePaths(movdir)

        for mov in movlist:
            basename = os.path.splitext(os.path.basename(mov))[0]
            outdir = os.path.join(stagedir, basename)
            mkdir_force(outdir)

            set_origins_and_spacing([fixed, mov])

            cmd = ['elastix',
                   '-f', fixed,
                   '-m', mov,
                   '-out', outdir,
                   '-p', elxparam_file,
                   '-threads', str(threads)]

            if fixed_mask:
                cmd.extend(['-fMask', fixed_mask])

            try:
                subprocess.check_output(cmd)
            except Exception:  # can't seem to log CalledProcessError
                logging.exception('registration falied:\n\n{}\n\n'.format(cmd))
                raise
            else:
                # Rename the registered output.
                elx_outfile = os.path.join(outdir, 'result.0.{}'.format(self.filetype))
                new_out_name = os.path.join(outdir, '{}.{}'.format(basename, self.filetype))
                shutil.move(elx_outfile, new_out_name)

                # Apply the transform to the mask so it can be used in the next stage
                tp_file = os.path.join(outdir, "TransformParameters.0.txt")
                #transform_mask(moving_mask, mask_dir_out, tp_file)

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
        except subprocess.CalledProcessError as e:
            sys.exit('### Transformix failed while transforming mask ###\nelastix command:{}'.format(cmd))
        else:
            os.remove(temp_tpfile)
            # Rename the transformed mask
            outname = os.path.join(mask_dir_out, 'result.nrrd')
            newname = os.path.join(mask_dir_out, os.path.basename(moving_mask))
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
        upper_extend = [d/2 for d in diffs]

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

        basename = os.path.splitext(os.path.basename(path))[0]
        padded_outname = os.path.join(outdir, '{}.{}'.format(basename, filetype))
        sitk.WriteImage(padded_vol, padded_outname, True)
    return 0


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
    spacing = (1.0, 1.0, 1.0)
    origin = (0.0, 0.0, 0.0)

    for vol in volpaths:
        im = sitk.ReadImage(vol)
        if im.GetSpacing() != spacing or im.GetOrigin() != origin:
            print 'setting orgin/spacing for {}'.format(vol)
            im.SetSpacing(spacing)
            im.SetOrigin(origin)
            sitk.WriteImage(im, vol)


def parse_yaml_config(configfile_or_dict):

    if not isinstance(configfile_or_dict, dict):  # Then it should be a yaml file
        try:
            config = yaml.load(open(configfile_or_dict, 'r'))
        except Exception as e:
            sys.exit("can't read the YAML config file - {}".format(e))
    else:
        config = configfile_or_dict
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

        sitk.WriteImage(average, out_path, True)   # Compressed=True
        #Check that it's been created
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
