#!/usr/bin/python

"""
..module:: The Medical Research Council registration pipeline
    :synopsis: Pipeline to create embryo average volumees.
.. moduleauthor:: Neil Horner
=================
Background
=================
The registration pipeline is for registereing volumes blah blah blah
There were some problems using masks with elastix 4.5. Version 4.7 worked fine

=================
Config file
=================
Currently works with config file v1

TODO:
    Override elastix parameters set in 'global elastix parameters' config section
    Make masks for each section go to output folder???
    Compress the masks for each stage. They will compress a lot
    Add linear interpolation to mask transformation
    Do not put deformation fields in  same directory as registration output

    If padding images, is the reference then updated to be the padded version?

"""

import os
import subprocess
import shutil
import sys
import argparse
import SimpleITK as sitk
import harwellimglib as hil
import yaml
import copy
import fnmatch
#import operator


def reg(configfile_or_dict):
    """This is the main function that is called by the GUI or from the command line.
    Reads in the config file, Creates directories, and initialises the registration process

    :param config: either a yaml file from the command line. Or a dict from the GUI
    :param: callback, function to return messages back to the GUI
    :return: nothing yet
    """
    config = parse_yaml_config(configfile_or_dict)

    config_report = validate_config(config)
    if config_report:
        sys.exit('\n\n'.join(config_report))

    filetype = config['filetype']

    # Get the variable from the config
    outdir = config['output_dir']
    if not os.path.isdir(outdir):
        sys.exit("The project directory: {} does not exit you have to create it".format(outdir))
    inputvols_dir = config['inputvolumes_dir']
    fixed_vol = config['fixed_volume']

    # Make the dir for the averages and the registration. If not set in config file, make default location
    if not 'average_dir' in config:
        average_dir = os.path.join(outdir, 'averages')
        config['average_dir'] = average_dir
    if not os.path.isdir(average_dir):
        os.mkdir(average_dir)

    # Pad out the moving dimensions of all the volumes to the same size. Do masks as well if required
    pad_dims = config.get('pad_dims')


    # TODO: Tidy this section up

    if pad_dims:
        paddedvols_dir = os.path.join(outdir, 'padded_inputs')
        mkdir_force(paddedvols_dir)
        pad_volumes(inputvols_dir, pad_dims, paddedvols_dir, config['filetype'])

    # Modify the config file to add the padded paths or specified vols to the first stage.
    first_stage_config = config['registration_stage_params'][0]
    fixed_vol_basname = os.path.basename(fixed_vol)

    if pad_dims:
        fixedvol_padded = os.path.splitext(fixed_vol_basname)[0]
        first_stage_config['fixed_vol'] = os.path.join(paddedvols_dir, '{}.{}'.format(fixedvol_padded, filetype))
        first_stage_config['movingvols_dir'] = paddedvols_dir
    else:
        first_stage_config['fixed_vol'] = fixed_vol
        first_stage_config['movingvols_dir'] = inputvols_dir

    do_registration(config)


def parse_yaml_config(configfile_or_dict):

    if not isinstance(configfile_or_dict, dict):  # Then it should be a yaml file
        try:
            config = yaml.load(open(configfile_or_dict, 'r'))
        except Exception as e:
            sys.exit("can't read the YAML config file - {}".format(e))
    else:
        config = configfile_or_dict
    return config


def validate_config(config):
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

    ### Check for correct paths
    if not os.path.isdir(config['inputvolumes_dir']):
        report.append("Input directory does not exit")

    if not os.path.isdir(config['output_dir']):
        report.append("Output directory does not exit")


    return report


def do_registration(config):
    """
    Run the registrations specified in the config file
    :param config: Config dictionary
    :return: 0 for success, or an error message
    """
    filetpye = config['filetype']
    outdir = config['output_dir']
    deformation_dir = os.path.join(outdir, 'deformation_fields')
    jacobians_dir = os.path.join(outdir, 'spatial_jacobians')

    delete_stages = []
    for i, reg_stage in enumerate(config['registration_stage_params']):

        #  Make the stage output dir
        stage_id = reg_stage['stage_id']
        stage_dir = os.path.join(outdir, stage_id)
        mkdir_force(stage_dir)

        if reg_stage.get('delete_stage_files'):
            delete_stages.append((i + 1, stage_dir))


        print("### Current registration step: {} ###".format(stage_id))

        # Make the elastix parameter file for this stage
        elxparam = make_elastix_parameter_file(config, reg_stage)
        elxparam_path = os.path.join(stage_dir, 'elastix_params_' + stage_id + '.txt')

        with open(elxparam_path, 'w') as fh:
            if elxparam:  # Not sure why I put this here
                fh.write(elxparam)

        # For rigid, miss out registering to self
        # Do a single stage registration
        elx_registration(elxparam_path,
                         reg_stage['fixed_vol'],
                         reg_stage['movingvols_dir'],
                         stage_dir,
                         config['threads'],
                         filetpye,
                         stage_id,
                         #reg_stage.get('mask_dir'],
                         #mask_dir_out,
                         reg_stage.get('fixed_mask')
                         )

        # Create a fixed mask from the transformed moving masks. ie. and avergae mask

        # Generate deformation fields, if required
        if reg_stage.get('generate_deformation_fields'):
            if not os.path.isdir(deformation_dir):
                os.mkdir(deformation_dir)
            if not os.path.exists(jacobians_dir):
                os.mkdir(jacobians_dir)
            generate_deformation_fields(stage_dir, deformation_dir, jacobians_dir, filetpye)

        # Make average
        average_path = os.path.join(config['average_dir'], '{0}.{1}'.format(stage_id, filetpye))
        make_average(stage_dir, average_path)

        for x in delete_stages:
            if x[0] == i:
                if os.path.isdir(x[1]):
                    shutil.rmtree(x[1])


        # Setup the fixed, moving and mask_dir for the next stage, if there is one
        if i + 1 < len(config['registration_stage_params']):
            next_reg_stage = config['registration_stage_params'][i + 1]

            if next_reg_stage.get('new_target'):
                next_reg_stage['fixed_vol'] = average_path  # The avergae from the previous step
            else:
                next_reg_stage['fixed_vol'] = reg_stage['fixed_vol'] # The same fixed volume as thje previous step

            if 'movingvols_dir' in next_reg_stage:  # Moving volumes are specified, don't just use previous stage
                stage_output_to_use = next_reg_stage['movingvols_dir']
                for stage in config['registration_stage_params']:
                    if stage['stage_id'] == stage_output_to_use:
                        next_reg_stage['movingvols_dir'] = os.path.join(outdir, stage_id)

            else:
                next_reg_stage['movingvols_dir'] = stage_dir  # The output dir of the previous registration

    print("### Registration pipeline finished ###")


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

def generate_deformation_fields(registration_dir, deformation_dir, jacobian_dir, filetype):
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
                deformation_out = os.path.join(deformation_dir, 'deformationField.{}'.format(filetype))
                jacobian_out = os.path.join(deformation_dir, 'spatialJacobian.{}'.format(filetype))
                deformation_renamed = os.path.join(deformation_dir, '{}_deformationFied.{}'.format(volid, filetype))
                jacobian_renamed = os.path.join(jacobian_dir, '{}_spatialJacobian.{}'.format(volid, filetype))
                shutil.move(deformation_out, deformation_renamed)
                shutil.move(jacobian_out, jacobian_renamed)


def make_elastix_parameter_file(config, stage_params):
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


def elx_registration(elxparam_file, fixed, movdir, stagedir, threads, filetype, stageid, fixed_mask):
    """
    Run elastix for one volume

    :return:
    """
    movlist= hil.GetFilePaths(movdir)

    for mov in movlist:
        basename = os.path.splitext(os.path.basename(mov))[0]
        outdir = os.path.join(stagedir, basename)
        mkdir_force(outdir)

        # Get the mask if it's require

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
        except subprocess.CalledProcessError as e:
            sys.exit('### Elastix failed ###\nCommand\n:{}'.format(e))
        else:
            # Rename the registered output.
            elx_outfile = os.path.join(outdir, 'result.0.{}'.format(filetype))
            new_out_name = os.path.join(outdir, '{}.{}'.format(basename, filetype))
            shutil.move(elx_outfile, new_out_name)

            # Apply the transform to the mask so it can be used in the next stage
            tp_file = os.path.join(outdir, "TransformParameters.0.txt")
            #transform_mask(moving_mask, mask_dir_out, tp_file)


def transform_mask(moving_mask, mask_dir_out, tp_file):
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


def pad_volumes(voldir, max_dims, outdir, filetype):
    """
    All volumes should be the same dimensions to aid in downstraem analysis and to prevent clipping when moving is
    larger than the fixed
    :return:
    """
    volpaths = hil.GetFilePaths(voldir)

    if len(volpaths) < 1:
        sys.exit("Can't find any volumes in {}".format(voldir))

    print 'Padding {} volumes/masks:'.format(str(len(volpaths)))

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


def logger():
    pass


def mkdir_force(dir_):
    if os.path.isdir(dir_):
        shutil.rmtree(dir_)
    os.mkdir(dir_)


if __name__ == "__main__":

    parser = argparse.ArgumentParser("The MRC Harwell image registration pipeline")
    parser.add_argument('-c', dest='config', help='Config file (YAML format)', required=True)
    args = parser.parse_args()

    #moving = hil.GetFilePaths(args.movdir)
    reg(args.config)
