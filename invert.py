#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""invert.py

This module inverts registrations performed with elastix

Example
-------

    $ invert.py -c invert.yaml

example config file:

    labelmap: padded_target/labelmap.nrrd
    voxel_size: 28
    stage_dirs:
    - deformable_to_8
    - deformable_to_128
    - affine
    - rigid

All paths are relative to the directory containing the config file


Notes
-----
The inversion will only work well for labelmaps as the final interpolation order is set to 0 to prevent interpolation
of lable map values and to keep them as the correct integers

Currently only inverts one elx_tform_params file per stage. Should be albe to do multple

Inversion can fail if the registration resolutions are set incorrectly.
For example, if the non-linear step has 6 resolutions and a a final BSpline grid spacing of 8, the largest grid size
will be 256. It seems that if this is larger than the input image dimensions, the inversion will fail.

"""

from os.path import join, splitext, abspath, basename
import os
import subprocess
import shutil
import sys
import yaml
import logging
from collections import defaultdict
import csv

import SimpleITK as sitk
import numpy as np

import harwellimglib as hil
import common

ELX_TRANSFORM_PREFIX = 'TransformParameters'
ELX_PARAM_PREFIX = 'elastix_params_'
INDV_REG_METADATA = 'reg_metadata.yaml'
FILE_FORMAT = '.nrrd'
LOG_FILE = 'inversion.log'
TRANSFORMIX_OUT = 'result.nrrd'
INVERSION_DIR_NAME = 'Inverted_transform_parameters'



def batch_invert_transform_parameters(config_file, outdir, threads=None):
    """
    stage_dirs: [list of stage directories to invert]

    Parameters
    ----------
    reg_dirs: list
        Ordered list of registration directories to apply inverse to. eg [deformable_to_8, deformable_to_128]
    outdir: str
        Absoulte path to output dir
    """
    with open(config_file, 'r') as yf:
        config = yaml.load(yf)

    config_dir = os.path.dirname(config_file)

    reg_dirs = list(reversed(config['reg_stages']))
    first_stage = join(config_dir, reg_dirs[0])
    volume_names = [basename(x) for x in hil.GetFilePaths(first_stage)]

    _mkdir_force(outdir)

    for i, vol_name in enumerate(volume_names):

        vol_id, vol_ext = splitext(vol_name)

        for r in reg_dirs:
            reg_dir = join(config_dir, r)

            #create output dir if not exists
            stage_out_dir = join(outdir, basename(reg_dir))
            if not os.path.exists(stage_out_dir):
                _mkdir_force(join(outdir, basename(reg_dir)))

            moving_dir = join(config_dir, reg_dir, vol_id)
            invert_param_dir = join(stage_out_dir, vol_id)


            stage_vol_files = os.listdir(moving_dir)  # All the files from the registration dir
            stage_files = os.listdir(reg_dir)         # The registration stage parent directory
            parameter_file = next(join(reg_dir, i) for i in stage_files if i.startswith(ELX_PARAM_PREFIX))
            transform_file = next(join(moving_dir, i) for i in stage_vol_files if i.startswith(ELX_TRANSFORM_PREFIX))

            if is_euler_stage(transform_file):
                continue  # We don't invert rigid/euler transforms using this method

            _mkdir_force(invert_param_dir)
            reg_metadata = yaml.load(open(join(moving_dir, INDV_REG_METADATA)))
            fixed_volume = join(moving_dir, reg_metadata['fixed_vol'])  # The original fixed volume used in the registration

            new_param = abspath(join(invert_param_dir, 'newParam.txt'))
            _modify_param_file(abspath(parameter_file), new_param)
            _invert_tform(fixed_volume, abspath(transform_file), new_param, invert_param_dir, threads)
            inverted_tform = abspath(join(invert_param_dir, 'TransformParameters.0.txt'))
            new_transform = abspath(join(invert_param_dir, 'inverted_transform.txt'))
            _modify_tform_file(inverted_tform, new_transform)




class BatchInvert(object):

    def __init__(self, config_path, threads=None):

        """
        Given a bunch of directories and the stages to invert to/from, do invert

        Parameters
        ----------
        config_path: str
            path to yaml config
        threads: str/ None
            number of threas to use. If None, use all available threads

        :return:
        """
        config_dir = os.path.dirname(config_path)
        logfile = os.path.join(config_dir, LOG_FILE)
        common.init_log(logfile, "Label inversion log")

        config = self.parse_yaml_config(config_path)
        self.voxel_size = config['voxel_size']
        self.organ_names = join(config_dir, config['organ_names'])

        self.volume_report = {}

        self.labelmap = join(config_dir, config['labelmap'])
        self.threads = threads
        self.out_dir = join(config_dir, 'inverted')
        common.mkdir_force(self.out_dir)

        self.volume_calc_path = join(self.out_dir, config['calculations_path'])

        self.registration_dirs = [join(config_dir, x) for x in config['stage_dirs']]
        self.elx_param_prefix = ELX_PARAM_PREFIX
        self.run()

    def parse_yaml_config(self, config_path):
        """
        Opens the yaml config file

        Parameters
        ----------
        config_path: str
            path to config file

        Returns
        -------
        dict:
            The config
        """

        try:
            config = yaml.load(open(config_path, 'r'))
        except Exception as e:
            sys.exit("can't read the YAML config file - {}".format(e))
        return config

    def run(self):
        """

        """

        volume_names = self.get_volume_names()
        organ_vol_calculations = defaultdict(dict)

        for i, vol_name in enumerate(volume_names):

            vol_id, vol_ext = splitext(vol_name)
            label_map = self.labelmap

            for reg_dir in self.registration_dirs:

                invert_stage_dir = join(self.out_dir, basename(reg_dir))
                if not os.path.isdir(invert_stage_dir):
                    os.mkdir(invert_stage_dir)

                invert_vol_dir = join(invert_stage_dir, vol_id)
                _mkdir_force(invert_vol_dir)

                moving_dir = join(reg_dir, vol_id)

                stage_vol_files = os.listdir(moving_dir)
                stage_files = os.listdir(reg_dir)
                parameter_file = next(join(reg_dir, i) for i in stage_files if i.startswith(self.elx_param_prefix))
                transform_file = next(join(moving_dir, i) for i in stage_vol_files if i.startswith(ELX_TRANSFORM_PREFIX))

                reg_metadata = yaml.load(open(join(moving_dir, INDV_REG_METADATA)))
                fixed = join(moving_dir, reg_metadata['fixed_vol'])

                inverted_label_map = process_volume(fixed, label_map, parameter_file, transform_file, invert_vol_dir, self.threads)
                if not inverted_label_map:
                    logging.warn('Inversion failed for {}'.format(reg_dir))

        calculate_organ_volumes(reg_dir, self.organ_names, self.volume_calc_path)

    def get_volume_names(self):
        """
        Volume names are obtained from the directory names from a single stage directory

        Returns
        -------
        List
            Volume names
        """
        first_stage = self.registration_dirs[0]
        volume_names = [basename(x) for x in hil.GetFilePaths(first_stage)]
        return volume_names


def process_volume(fixed_volume, invertable_volume, parameter_file, transform_file, outdir, threads):
    """The main function for doing the inversion

    Parameters
    ----------
    fixed_volume : str
        fixed volume path
    invertable_volume : str
        volume to invert
    parameter_file: str
        path to elastix parameter file that was used in the original registration
    transform_file: str
        path to the elastix output transform paramter file
    outdir: str
        path to place output from this function
    threads: str
        number of threads to use for elastix and transformix

    Returns
    -------
    str/bool
        path to new img if succesful else False

    Raises
    ------
    subprocess.CalledProcessError
        If the elastix or transformix subprocess call fails
    """

    _mkdir_force(outdir)

    new_param = abspath(join(outdir, 'newParam.txt'))
    _modify_param_file(abspath(parameter_file), new_param)
    _invert_tform(fixed_volume, abspath(transform_file), new_param, outdir, threads)
    inverted_tform = abspath(join(outdir, 'TransformParameters.0.txt'))
    new_transform = abspath(join(outdir, 'inverted_transform.txt'))
    _modify_tform_file(inverted_tform, new_transform)

    specimen_basename = os.path.basename(outdir)
    new_img_path = os.path.join(outdir, 'seg_' + specimen_basename + FILE_FORMAT)
    tform_success = _invert_image(invertable_volume, new_transform, outdir, new_img_path, threads)
    if tform_success:

        return new_img_path
    else:
        return False


def _modify_param_file(elx_param_file, newfile_name):
    """
    Modifies the elastix input parameter file that was used in the original transformation.
    Adds DisplacementMagnitudePenalty (which is needed for inverting)
    Turns off writing the image results at the end as we only need an inveterted output file

    Parameters
    ----------
    elx_param_file: str
        path to elastix input parameter file
    newfile_name: str
        path to save modified parameter file to

    """
    with open(elx_param_file) as old, open(newfile_name, "w") as new:

        for line in old:
            if line.startswith("(Metric"):
                line = '(Metric "DisplacementMagnitudePenalty")\n'
            if line.startswith('(WriteResultImage'):
                line = '(WriteResultImage "false")\n'
            if line.startswith('(FinalBSplineInterpolationOrder'):
                line = '(FinalBSplineInterpolationOrder  0)\n'
            if line.startswith('(RigidityPenaltyWeight'):  # a test just for rigidity penalty
                line = ''
            new.write(line)


def _invert_tform(fixed, tform_file, param, outdir, threads):
    """
    Invert the transform and get a new transform file
    """
    cmd = ['elastix',
           '-t0', tform_file,
           '-p', param,
           '-f', fixed,
           '-m', fixed,
           '-out', outdir
           ]
    if threads:
        cmd.extend(['-threads', threads])
    try:
        subprocess.check_output(cmd)
    except Exception as e:
        print 'elastix failed: {}'.format(e)
        logging.error('Inverting transform file failed. cmd: {}:\nmessage'.format(cmd), exc_info=True)



def _modify_tform_file(inverted_tform, newfile_name):
    """
    Remove "NoInitialTransform" from the output transform parameter file
    Set output image format to unsigned char
    args:
        outdir: directory where the inverted transform file is
        newfile_mame: where to save modified transform file
    """

    new_tform_param_fh = open(newfile_name, "w")
    tform_param_fh = open(inverted_tform, "r")
    for line in tform_param_fh:
        if line.startswith('(InitialTransformParametersFileName'):
            line = '(InitialTransformParametersFileName "NoInitialTransform")\n'
        if line.startswith('(ResultImagePixelType'):
            line = '(ResultImagePixelType "unsigned char")\n'
        new_tform_param_fh.write(line)
    new_tform_param_fh.close()
    tform_param_fh.close()


def _invert_image(vol, tform, outdir, rename_output, threads):
    """
    Using the iverted elastix transform paramter file, invert a volume with transformix

    Parameters
    ----------
    vol: str
        path to volume to invert
    tform: str
        path to elastix transform parameter file
    outdir: str
        path to save transformix output
    rename_output: str
        rename the transformed volume to this
    threads: str/None
        number of threads for transformix to use. if None, use all available cpus
    Returns
    -------
    str/bool
        path to new img if succesful else False
    """
    cmd = ['transformix',
           '-in', vol,
           '-tp', tform,
           '-out', outdir
           ]

    if threads:
        cmd.extend(['-threads', threads])
    try:
        subprocess.check_output(cmd)
    except Exception as e:
        print 'transformix failed'
        logging.error('transformix failed with this command: {}\nerror message:'.format(cmd), exc_info=True)
        return False
    try:
        old_img = os.path.join(outdir, TRANSFORMIX_OUT)
        os.rename(old_img, rename_output)
    except OSError:
        return False
    else:
        return True


def calculate_organ_volumes(labels_dir, label_names, outfile, voxel_size=28.0):
    """
    Calculate organ volumes and save results in a csv file

    Parameters
    ----------
    labels_dir: str
        path to dir containing label maps. Volumes can exist in subdirectories.
    label_names: str
        path to csv file containing label names:
        eg:
            liver
            thyroid
            left kidney
    outfile: str
        path to save calculated volumes as csv file
    voxel_size: float
        size of voxels
    """
    header = [' ']
    with open(label_names, 'r') as headerfile:
        reader = csv.reader(headerfile)
        for row in reader:
            header.append(row[0])

    with open(outfile, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(header)
        labelmaps = hil.GetFilePaths(labels_dir)
        for label_path in labelmaps:
            full_path = join(labels_dir, label_path)
            labelmap = sitk.ReadImage(full_path)
            label_array = sitk.GetArrayFromImage(labelmap)
            hist = np.histogram(label_array, bins=label_array.max() + 1)[0]
            row = list(hist)
            #Remove first value. as this is the background and replace with vol name
            row[0] = basename(label_path)
            # Need to normalise to voxel size
            writer.writerow(row)

def is_euler_stage(tform_param):
    """
    Return True if the registration used to create this param file was a Euler transform. Can't currently invert
    Euler transforms with this method, and is usually not required
    :param tform_param:
    :return:
    """
    with open(tform_param, 'r') as fh:
        line = fh.readline()
        if 'EulerTransform' in line:
            return True
        else:
            return False

def _mkdir_force(dir_):
    if os.path.isdir(dir_):
        shutil.rmtree(dir_)
    os.mkdir(dir_)


if __name__ == '__main__':
    import argparse

    if sys.argv[1] == 'labels':
        parser = argparse.ArgumentParser("calculate organ volumes")
        parser.add_argument('-l', '--label_dir', dest='labels', help='directory with labels', required=True)
        parser.add_argument('-n', '--label_names', dest='label_names', help='csv files with label names', required=True)
        parser.add_argument('-o', '--out', dest='outfile', help='path to csv file', required=True)
        #parser.add_argument('-n', '--label_names', dest='label_names', help='csv files with label names', required=True)
        args, _ = parser.parse_known_args()
        calculate_organ_volumes(args.labels, args.label_names, args.outfile)

    elif sys.argv[1] == 'invert':
        parser = argparse.ArgumentParser("invert elastix registrations and calculate organ volumes")
        parser.add_argument('-c', '--config', dest='config', help='yaml config file', required=True)
        parser.add_argument('-t', '--threads', dest='threads', type=str, help='number of threads to use', required=False)

        args, _ = parser.parse_known_args()

    elif sys.argv[1] == 'invert_reg':
        parser = argparse.ArgumentParser("invert elastix registrations and calculate organ volumes")
        parser.add_argument('-c', '--config',  dest='config', help='Config file with list of registration dirs', required=True)
        parser.add_argument('-o', '--out',  dest='outdir', help='where to put the output', required=True)

        parser.add_argument('-t', '--threads', dest='threads', type=str, help='number of threads to use', required=False)

        args, _ = parser.parse_known_args()
        batch_invert_transform_parameters(args.config, args.outdir)


        #
        # if not args.threads:
        #     args.threads = None
        #
        # BatchInvert(args.config, args.threads)
