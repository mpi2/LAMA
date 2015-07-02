#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""invert.py

This module inverts registrations performed with elastix

Example
-------

    $ invert.py -f fixed.nrrd, -m moving.nrrd -p original_elx-file -t elx-tform-params

Notes
-----
The inversion will only work well for labelmaps as the final interpolation order is set to 0 to keep the
label map values intact
Currently only inverts one elx_tform_params file per stage. Should be albe to do multple

TODO: Currently in batch mode, the inversions end up in wrong dirs. 128 in 8, affine in 128 etc

"""

from os.path import join, splitext, abspath, basename
import os
import subprocess
import shutil
import yaml
import logging
import common

ELX_TRANSFORM_PREFIX = 'TransformParameters'
INDV_REG_METADATA = 'reg_metadata.yaml'
FILE_FORMAT = '.nrrd'
LOG_FILE = 'inversion.log'
TRANSFORMIX_OUT = 'result.nrrd'


class BatchInvert(object):
    def __init__(self, label_map, stage_dirs, volume_names, out_dir, elx_param_prefix, threads):
        """
        Given a bunch of directories and the stages to invert to/from, do invert

        Parameters
        ----------
        label_map: str
            path to label map or other data to be inverted
        out_dir: str
            path to save inversions to
        stage_dirs: iterable
            sequence of registration directories. In the order of initial registration
        volume_names: list
            basenames of volumes
        elx_param_prefix: str
            the prefix added to elastix paramter files
        threads: str
            number of threads to use for tranformations

        :return:
        """
        logfile = os.path.join(out_dir, LOG_FILE)
        common.init_log(logfile, "Label inversion log")

        self.labelmap = label_map
        self.threads = threads
        self.out_dir = out_dir
        self.volume_names = volume_names
        self.registration_dirs = stage_dirs
        self.elx_param_prefix = elx_param_prefix
        self.run()

    def run(self):

        for vol_name in self.volume_names:
            vol_id, vol_ext = splitext(vol_name)
            label_map = self.labelmap
            for reg_dir in reversed(self.registration_dirs):
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

                label_map = process_volume(fixed, label_map, parameter_file, transform_file, invert_vol_dir, self.threads)
                if not label_map:
                    logging.warn('Inversion failed for {}'.format(reg_dir))


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

    """
    with open(elx_param_file) as old, open(newfile_name, "w") as new:

        for line in old:
            if line.startswith("(Metric"):
                line = '(Metric "DisplacementMagnitudePenalty")\n'
            if line.startswith('(WriteResultImage'):
                line = '(WriteResultImage "false")\n'
            if line.startswith('(FinalBSplineInterpolationOrder'):
                line = '(FinalBSplineInterpolationOrder  0)\n'
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
           '-out', outdir,
           '-threads', threads
           ]
    try:
        subprocess.check_output(cmd)
    except Exception as e:
        print 'elastix failed: {}'.format(e)
        logging.error('Inverting transform file failed. cmd: {}:\nmessage {}'.format(cmd, e))



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
    Parameters
    ----------

    Returns
    -------
    str/bool
        path to new img if succesful else False
    """
    cmd = ['transformix',
           '-in', vol,
           '-tp', tform,
           '-out', outdir,
           '-threads', threads
           ]
    try:
        subprocess.check_output(cmd)
    except Exception as e:
        print 'transformix failed'
        logging.error('transformix failed with this command: {}\nerror message:{}'.format(cmd, e))
        return False
    try:
        old_img = os.path.join(outdir, TRANSFORMIX_OUT)
        os.rename(old_img, rename_output)
    except OSError:
        return False
    else:
        return True




def _mkdir_force(dir_):
    if os.path.isdir(dir_):
        shutil.rmtree(dir_)
    os.mkdir(dir_)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser("invert elastix registrations")
    parser.add_argument('-f', '--fixed', dest='fixed', help='fixed volume ', required=True)
    parser.add_argument('-m', '--moving', dest='moving', help='moving volume', required=True)
    parser.add_argument('-p', '--param', dest='param', help='original elastix parameter file used for registration', required=True)
    parser.add_argument('-t', '--tform', dest='tform', help='elastix transform parameter output filr', required=True)
    parser.add_argument('-o', '--out', dest='outdir', help='', required=True)
    parser.add_argument('-th', '--threads', dest='threads', help='', type=str, default='1')

    args = parser.parse_args()

    process_volume(args.fixed, args.moving, args.param, args.tform, args.outdir, args.threads)
