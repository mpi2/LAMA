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
ELX_INVERTED_POINTS_NAME = 'outputpoints.vtk'
INDV_REG_METADATA = 'reg_metadata.yaml'
FILE_FORMAT = '.nrrd'
LOG_FILE = 'inversion.log'
TRANSFORMIX_OUT = 'result.nrrd'
INVERSION_DIR_NAME = 'Inverted_transform_parameters'
INVERTED_TRANSFORM_NAME = 'inverted_transform.txt'
VOLUME_CALCULATIONS_FILENAME = "organvolumes.csv"


def batch_invert_transform_parameters(config_file, invert_config_file, outdir, threads=None):
    """
    Invert registrations creating new TransformParameter files that can then be used by transformix to invert labelmaps
     etc

    Parameters
    ----------
    config_file: str
        path to config file. Can be manually created or is generated by mrch_regpipeline.py. Contains paths to registered
        image directories. In order of whuch they were registered

        Example config:

           registration_stage_params:
                - Affine
                - deformable_128
                - deformable_16
                - deformable_8

    outdir: str
        Absoulte path to output dir

    invert_config_file: str
        path to output inverion config to
    """

    with open(config_file, 'r') as yf:
        config = yaml.load(yf)
        config_dir = os.path.dirname(config_file)

    reg_dirs = get_reg_dirs(config, config_dir)

    first_stage = join(config_dir, reg_dirs[0])
    volume_names = [basename(x) for x in hil.GetFilePaths(first_stage)]

    _mkdir_force(outdir)
    stages_to_invert = defaultdict(list)

    for i, vol_name in enumerate(volume_names):

        vol_id, vol_ext = splitext(vol_name)

        for r in reg_dirs:
            reg_dir = join(config_dir, r)

            #create output dir if not exists
            stage_out_dir = join(outdir, basename(reg_dir))

            moving_dir = join(config_dir, reg_dir, vol_id)
            invert_param_dir = join(stage_out_dir, vol_id)

            stage_vol_files = os.listdir(moving_dir)  # All the files from the registration dir
            stage_files = os.listdir(reg_dir)         # The registration stage parent directory
            parameter_file = next(join(reg_dir, i) for i in stage_files if i.startswith(ELX_PARAM_PREFIX))
            transform_file = next(join(moving_dir, i) for i in stage_vol_files if i.startswith(ELX_TRANSFORM_PREFIX))

            if is_euler_stage(transform_file):
                continue  # We don't invert rigid/euler transforms using this method

            if not os.path.exists(stage_out_dir):
                _mkdir_force(join(outdir, basename(reg_dir)))

            rel_inversion_path = os.path.basename(r)
            if rel_inversion_path not in stages_to_invert['inversion_order']:
                stages_to_invert['inversion_order'].insert(0, rel_inversion_path)
            _mkdir_force(invert_param_dir)
            reg_metadata = yaml.load(open(join(moving_dir, INDV_REG_METADATA)))
            fixed_volume = join(moving_dir, reg_metadata['fixed_vol'])  # The original fixed volume used in the registration

            new_param = abspath(join(invert_param_dir, 'newParam.txt'))
            _modify_param_file(abspath(parameter_file), new_param)
            _invert_tform(fixed_volume, abspath(transform_file), new_param, invert_param_dir, threads)
            inverted_tform = abspath(join(invert_param_dir, 'TransformParameters.0.txt'))
            new_transform = abspath(join(invert_param_dir, 'inverted_transform.txt'))
            _modify_tform_file(inverted_tform, new_transform)

    # Create a yaml config file so that inversions can be run seperatley
    with open(invert_config_file, 'w') as yf:
        yf.write(yaml.dump(dict(stages_to_invert), default_flow_style=False))


def get_reg_dirs(config, config_dir):
    """

    """

    reg_stages = []
    for i, reg_stage in enumerate(config['registration_stage_params']):
        stage_id = reg_stage['stage_id']
        stage_dir = join(config_dir, config['output_dir'], stage_id)
        reg_stages.append(stage_dir)
    return reg_stages


class Invert(object):
    def __init__(self, config_path, invertable_volume, outdir, threads=None):
        """
        Inverts a bunch of volumes/label maps. A yaml config file specifies the order of inverted transform parameters
        to use. This config file should be in the root of the directory containing these inverted tform dirs.

        Also need to input a directory containing volumes/label maps etc to invert. These need to be in directories
        named with the same name as the corresponding inverted tform file directories

        Parameters
        ----------
        config_path: str
            path to yaml config containing the oder of the inverted directories to use
        threads: str/ None
            number of threas to use. If None, use all available threads
        invert_single: bool
            if true only invert using invert paraneters with the sdame basename. EG for when inverting stats

        :return:
        """
        with open(config_path, 'r') as yf:
            self.config = yaml.load(yf)

        self.config_dir = os.path.dirname(config_path)
        
        self.invertable_object = invertable_volume
        self.threads = threads
        self.out_dir = outdir
        common.mkdir_if_not_exists(self.out_dir)

        self.inverted_tform_stage_dirs = self.get_inversion_dirs()

        self.invertable_basename = os.path.splitext(os.path.basename(invertable_volume))[0]

        self.elx_param_prefix = ELX_PARAM_PREFIX

        self.run()

    def get_inversion_dirs(self):

        dirs = []
        for dir_name in self.config['inversion_order']:
            dir_path = join(self.config_dir, dir_name)
            dirs.append(dir_path)

        return dirs

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

        for inversion_stage in self.inverted_tform_stage_dirs:
            invert_stage_out = join(self.out_dir, basename(inversion_stage))
            if not os.path.isdir(invert_stage_out):
                common.mkdir_if_not_exists(invert_stage_out)

            inv_tform_dir = join(inversion_stage, self.invertable_basename)

            transform_file = join(inv_tform_dir, INVERTED_TRANSFORM_NAME)
            if os.path.exists(transform_file):
                print "Error. Cannot find inversion transform file: {}".format(transform_file)
                return
            invert_vol_out_dir = join(invert_stage_out, self.invertable_basename)
            common.mkdir_if_not_exists(invert_vol_out_dir)

            invertable = self._invert(invertable, transform_file, invert_vol_out_dir, self.threads)

    def _invert(self):
        raise NotImplementedError


class InvertLabelMap(Invert):

    def __init__(self, *args, **kwargs):
        super(InvertLabelMap, self).__init__(*args, **kwargs)
        a = args
        k = kwargs

    def run(self):
        """
        Calls the parent run function to invert the labels.
        Then optionally calculates organ volumes for the final inverted labels
        """
        super(InvertLabelMap, self).run()

        if self.do_organ_vol_calcs:
            # get the last invert stage. And calculate organ volumes for this
            last_invert_stage = self.inverted_tform_stage_dirs[-1]
            invert_stage_out = join(self.out_dir, basename(last_invert_stage))
            organ_size_out = join(self.out_dir, VOLUME_CALCULATIONS_FILENAME)
            calculate_organ_volumes(invert_stage_out, self.organ_names, organ_size_out)

    def _invert(self, labelmap, tform, outdir, rename_output, threads=None):
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
        lm_basename = os.path.splitext(os.path.basename(labelmap))[0]
        new_img_path = join(outdir, 'seg_' + lm_basename + FILE_FORMAT)

        cmd = [
            'transformix',
            '-in', labelmap,
            '-tp', tform,
            '-out', outdir
        ]

        if threads:
            cmd.extend(['-threads', threads])
        try:
            subprocess.check_output(cmd)
        except Exception as e:
            print 'transformix failed inverting labelmap: {}'.format(labelmap)
            #logging.error('transformix failed with this command: {}\nerror message:'.format(cmd), exc_info=True)
            return False
        try:
            old_img = os.path.join(outdir, TRANSFORMIX_OUT)
            os.rename(old_img, new_img_path)
        except OSError:

            return old_img
        else:
            return new_img_path


class InvertMeshes(Invert):

    def __init__(self, *args, **kwargs):
        super(InvertMeshes, self).__init__(*args, **kwargs)

    def _invert(self, mesh, tform, outdir, threads=None):
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
        m_basename = os.path.splitext(os.path.basename(mesh))[0]
        new_vtk_path = join(outdir, m_basename + '.vtk')

        cmd = [
            'transformix',
            '-def', mesh,
            '-tp', tform,
            '-out', outdir
        ]

        if threads:
            cmd.extend(['-threads', threads])
        try:
            subprocess.check_output(cmd)
        except Exception as e:
            print 'transformix failed inverting mesh: {}'.format(mesh)
            #logging.error('transformix failed with this command: {}\nerror message:'.format(cmd), exc_info=True)
            return False
        try:
            old_vtk = os.path.join(outdir, ELX_INVERTED_POINTS_NAME)
            os.rename(old_vtk, new_vtk_path)
        except OSError:

            return old_vtk
        else:
            return new_vtk_path

class InvertVol(Invert):
    def __del__(self, *args):
        super(InvertVol, self).__init__(*args)

    def _invert(self, volume, tform, outdir, threads=None):
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

        lm_basename = os.path.splitext(os.path.basename(volume))[0]
        new_img_path = join(outdir, 'seg_' + lm_basename + FILE_FORMAT)

        cmd = [
            'transformix',
            '-in', volume,
            '-tp', tform,
            '-out', outdir
        ]

        if threads:
            cmd.extend(['-threads', threads])
        try:
            subprocess.check_output(cmd)
        except Exception as e:
            print 'transformix failed inverting labelmap: {}'.format(volume)
            #logging.error('transformix failed with this command: {}\nerror message:'.format(cmd), exc_info=True)
            return False
        try:
            old_img = os.path.join(outdir, TRANSFORMIX_OUT)
            os.rename(old_img, new_img_path)
        except OSError:

            return old_img
        else:
            return new_img_path



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
            if line.startswith("(Metric "):
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
        #logging.error('Inverting transform file failed. cmd: {}:\nmessage'.format(cmd), exc_info=True)


def _modify_tform_file(elx_tform_file, newfile_name):
    """
    Remove "NoInitialTransform" from the output transform parameter file
    Set output image format to unsigned char. Writes out a modified elastix transform parameter file
    that can be used for inverting volumes

    Parameters
    ----------
    elx_tform_file: str
        path to elastix transform file
    newfile_mame: str
        path to save modified transform file
    """

    new_tform_param_fh = open(newfile_name, "w")
    tform_param_fh = open(elx_tform_file, "r")
    for line in tform_param_fh:
        if line.startswith('(InitialTransformParametersFileName'):
            line = '(InitialTransformParametersFileName "NoInitialTransform")\n'
        if line.startswith('(ResultImagePixelType'):
            line = '(ResultImagePixelType "unsigned char")\n'
        new_tform_param_fh.write(line)
    new_tform_param_fh.close()
    tform_param_fh.close()


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
    if label_names:
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
            # Get the name of the volume
            volname = os.path.split(os.path.split(full_path)[0])[1]
            # Remove first value. as this is the background and replace with vol name
            row[0] = basename(volname)
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

    elif sys.argv[1] == 'invert_labels':
        parser = argparse.ArgumentParser("invert elastix registrations and calculate organ volumes")
        parser.add_argument('-c', '--config', dest='config', help='yaml config file', required=True)
        parser.add_argument('-i', '--invertable', dest='invertable', help='volume to invert', required=True)
        parser.add_argument('-o', '--outdir', dest='outdir', help='output dir', required=True)
        parser.add_argument('-t', '--threads', dest='threads', type=str, help='number of threads to use', required=False)

        args, _ = parser.parse_known_args()
        InvertLabelMap(args.config, args.invertable, args.outdir)

    elif sys.argv[1] == 'invert_reg':
        parser = argparse.ArgumentParser("invert elastix registrations to create elastix inversion parameter files")
        parser.add_argument('-c', '--config',  dest='config', help='Config file with list of registration dirs', required=True)
        parser.add_argument('-o', '--out',  dest='outdir', help='where to put the output', required=True)
        parser.add_argument('-t', '--threads', dest='threads', type=str, help='number of threads to use', required=False)
        args, _ = parser.parse_known_args()
        batch_invert_transform_parameters(args.config, args.outdir)

    elif sys.argv[1] == 'invert_vol':
        parser = argparse.ArgumentParser("invert elastix registrations and calculate organ volumes")
        parser.add_argument('-c', '--config', dest='config', help='yaml config file', required=True)
        parser.add_argument('-i', '--invertable', dest='invertable', help='volume to invert', required=True)
        parser.add_argument('-o', '--outdir', dest='outdir', help='output dir', required=True)
        parser.add_argument('-t', '--threads', dest='threads', type=str, help='number of threads to use', required=False)

        args, _ = parser.parse_known_args()
        Invert(args.config, args.invertable, args.outdir)


        #
        # if not args.threads:
        #     args.threads = None
        #
        # BatchInvert(args.config, args.threads)
