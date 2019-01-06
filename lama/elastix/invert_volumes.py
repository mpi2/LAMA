#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""invert_volumes.py

This module inverts registrations performed with elastix

Example
-------

    $ invert_values.py -c invert.yaml

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
from pathlib import Path
import os
import subprocess
from os.path import join, abspath, isfile
from typing import Union

from logzero import logger as logging
import yaml

from lama import common
from lama.registration_pipeline.validate_config import LamaConfig

from lama.elastix.invert_transforms import \
ELX_TRANSFORM_PREFIX, ELX_PARAM_PREFIX, ELX_INVERTED_POINTS_NAME, TRANSFORMIX_OUT_NAME, \
LABEL_INVERTED_TRANFORM, IMAGE_INVERTED_TRANSFORM

common.add_elastix_env()



class Invert(object):
    def __init__(self, config_path, invertable, outdir, threads=None, noclobber=False):
        """
        Inverts a series of volumes. A yaml config file specifies the order of inverted transform parameters
        to use. This config file should be in the root of the directory containing these inverted tform dirs.

        Also need to input a directory containing volumes/label maps etc to invert. These need to be in directories
        named with the same name as the corresponding inverted tform file directories

        Parameters
        ----------
        config_path: str
            path to yaml config containing the oder of the inverted directories to use
        threads: str/ None
            number of threas to use. If None, use all available threads
        invertable_volume: str
            path to object to invert
        invertable: str
            dir or path. If dir, invert all objects within the subdirectories.
                If path to object (eg. labelmap) invert that instead
        noclobber: bool
            if True do not overwrite already inverted labels
        :return:
        """

        self.noclobber = noclobber

        with open(config_path, 'r') as yf:
            self.config = yaml.load(yf)

        self.invertables = invertable
        self.config_dir = os.path.dirname(config_path)  # The dir containing the inverted elx param files

        self.threads = threads
        self.out_dir = outdir
        common.mkdir_if_not_exists(self.out_dir)

        self.inverted_tform_stage_dirs = self.get_inversion_dirs()
        self.forward_tform_stage_dirs = self.get_forward_tranforms()

        self.elx_param_prefix = ELX_PARAM_PREFIX
        self.invert_transform_name = None  # Set in subclasses
        self.last_invert_dir = None

    def get_inversion_dirs(self):

        dirs = []
        for dir_name in self.config['inversion_order']:
            dir_path = join(self.config_dir, dir_name)
            dirs.append(dir_path)
        return dirs

    def get_forward_tranforms(self):
        dirs = []
        reg_dir = self.config.get('registration_directory')
        for dir_name in self.config['inversion_order']:
            dir_path = join(self.config_dir, reg_dir, dir_name)
            dirs.append(dir_path)
        return dirs


    @staticmethod
    def parse_yaml_config(config_path):
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
        done_file = abspath(join(self.out_dir, 'invert.done'))

        if not isfile(done_file):
            open(done_file, 'a').close()

        # lock_path = f'{done_file}.lock'
        #
        # lock = FileLock(lock_path, timeout=20)

        inverting_names = os.listdir(self.inverted_tform_stage_dirs[0])

        for i, vol_name in enumerate(inverting_names):

            # Get a lock on the progress log, check if the vol(specimen) has been done
            # with lock:
            #     with open(done_file, 'r+') as fh:
            #
            #         done = [x.strip() for x in fh.readlines()]
            #         if vol_name in done:
            #             print(f'skipping {vol_name}')
            #             continue
            #         else:
            #             print(f'inverting {vol_name}')
            #             fh.write(f'{vol_name}\n')
            #             fh.flush()

            invertable = self.invertables

            for inversion_stage, forward_stage in zip(self.inverted_tform_stage_dirs, self.forward_tform_stage_dirs):
                invert_stage_out = join(self.out_dir, basename(inversion_stage))
                if not os.path.isdir(invert_stage_out):
                    common.mkdir_if_not_exists(invert_stage_out)

                if self.type == 'forward':  # temp bodge for mesh inversion problem
                    inv_tform_dir = join(forward_stage, vol_name)
                    transform_file = join(inv_tform_dir, self.invert_transform_name)
                else:
                    inv_tform_dir = join(inversion_stage, vol_name)
                    transform_file = join(inv_tform_dir, self.invert_transform_name)

                invert_vol_out_dir = join(invert_stage_out, vol_name)

                common.mkdir_if_not_exists(invert_vol_out_dir)

                logging.info('inverting {}'.format(transform_file))

                invertable = self._invert(invertable, transform_file, invert_vol_out_dir, self.threads)

                if not invertable: # If inversion failed or there is nocobber, will get None
                    continue # Move on to next volume to invert
        # lock.release()

        self.last_invert_dir = invert_stage_out

    def _invert(self):
        raise NotImplementedError


class InvertLabelMap(Invert):

    def __init__(self, *args, **kwargs):
        super(InvertLabelMap, self).__init__(*args, **kwargs)
        self.invert_transform_name = LABEL_INVERTED_TRANFORM
        self.type = 'normal'

    def run(self):
        """
        Calls the parent run function to invert the labels.
        Then optionally calculates organ volumes for the final inverted labels
        """
        super(InvertLabelMap, self).run()

    def _invert(self, labelmap, tform, outdir, threads=None):
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
        #lm_basename = os.path.splitext(os.path.basename(labelmap))[0]
        if not common.test_installation('transformix'):
            raise OSError('Cannot find transformix. Is it installed?')

        old_img = os.path.join(outdir, TRANSFORMIX_OUT)                # where thetransformix-inverted labelmap will be

        path, base = os.path.split(os.path.normpath(outdir))
        new_output_name = os.path.join(outdir, '{}.nrrd'.format(base)) # Renamed transformix-inverted labelmap

        # if self.noclobber and isfile(new_output_name):
        # Maybe need to do two types of noclobber
        # 1: where if the folder does not exist, do not do it
        # 2: where the folder exists but the final output file does not exist
        #     return None


        cmd = [
            'transformix',
            '-in', labelmap,
            '-tp', tform,
            '-out', outdir
        ]

        if threads:
            cmd.extend(['-threads', str(threads)])
        try:
            subprocess.check_output(cmd)
        except Exception as e:
            logging.exception('{}\ntransformix failed inverting labelmap: {}'.format(e, labelmap))
            # sys.exit()
            logging.error('transformix failed with this command: {}\nerror message:'.format(cmd))

        try:
            shutil.move(old_img, new_output_name)
        except IOError as e:
            print()
            'could not rename {}'.format(old_img)
            return old_img
        else:
            return new_output_name


class InvertStats(InvertLabelMap):
    """
    This class behaves almost the same as InvertLabelMap in that it inverts a single image file back onto multiple
    inputs. It just uses a different elastix parameters
    """
    def __init__(self, *args, **kwargs):
        super(InvertStats, self).__init__(*args, **kwargs)
        self.invert_transform_name = IMAGE_INVERTED_TRANSFORM
        self.type = 'normal'


class InvertMeshes(Invert):

    def __init__(self, config_path, invertable, outdir, threads=None):
        super(InvertMeshes, self).__init__(config_path, invertable, outdir, threads)
        self.invert_transform_name = ELX_TRANSFORM_PREFIX
        self.type = 'forward'

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
        common.test_installation('transformix')
        m_basename = os.path.splitext(os.path.basename(mesh))[0]
        new_vtk_path = join(outdir, m_basename + '.vtk')

        cmd = [
            'transformix',
            '-def', mesh,
            '-tp', tform,
            '-out', outdir
        ]

        if threads:
            cmd.extend(['-threads', str(threads)])
        try:
            subprocess.check_output(cmd)
        except Exception as e:
            print('transformix failed inverting mesh: {}'.format(mesh))
            logging.error('transformix failed with this command: {}\nerror message:'.format(cmd), exc_info=True)
            print(e)
            sys.exit(1)
        try:
            # rename the inverted points form this stage
            old_vtk = os.path.join(outdir, ELX_INVERTED_POINTS_NAME)
            os.rename(old_vtk, new_vtk_path)
        except OSError:

            raise
        else:
            return new_vtk_path


# class InvertRoi(InvertLabelMap):
#     def __init__(self, config_path, invertable, outdir, vol_info, voxel_size, threads=None):
#         super(InvertRoi, self).__init__(config_path, invertable, outdir, threads)
#         self.invert_transform_name = LABEL_INVERTED_TRANFORM
#         self.vol_info = vol_info
#         self.voxel_size = voxel_size
#
#     def run(self):
#         super(InvertRoi, self).run()
#         # At this point we have a bunch of rois inverted onto the padded inputs
#         # We need to adjust the rois to account for the padding
#         out = join(self.out_dir, 'Extracted_roi')
#         unpad_roi(self.vol_info, self.last_invert_dir, self.voxel_size, out)


class InvertSingleVol(Invert):
    """
    Invert volumes using the elastix inverted transform parameters.
    This class is used for inverting statitistics overlays
    """
    def __init__(self, *args, **kwargs):
        super(InvertSingleVol, self).__init__(*args, **kwargs)
        self.invert_transform_name = IMAGE_INVERTED_TRANSFORM

    def run(self, prefix=None):
        """
        Parameters
        ----------
        prefix: str
            A prefix that is added to the stats volumes. To locate correct transform inversion files, look for files
            with this prefix missing
        """

    #     inverting_names = os.listdir(self.inverted_tform_stage_dirs[        # 0])

    #     for i, vol_name in enumerate(inverting_names):
    #         if self.batch_invert:
    #             invertable = self.invertables
    #         else:
    #             invertable = self.invertables[vol_name]

        volname, ext = splitext(basename(self.invertables))
        if prefix and volname.startswith(prefix):
            original_vol_name = volname[len(prefix):]  # remove the prfix to get the original vol name to find the tf
        else:
            original_vol_name = volname
        invertable = self.invertables

        for inversion_stage in self.inverted_tform_stage_dirs:
            invert_stage_out = join(self.out_dir, basename(inversion_stage))
            common.mkdir_if_not_exists(invert_stage_out)

            inv_tform_dir = join(inversion_stage, original_vol_name)

            transform_file = join(inv_tform_dir, IMAGE_INVERTED_TRANSFORM)
            invert_vol_out_dir = join(invert_stage_out, volname)
            common.mkdir_if_not_exists(invert_vol_out_dir)

            invertable = self._invert(invertable, transform_file, invert_vol_out_dir, self.threads)

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
        common.test_installation('transformix')
        lm_basename = os.path.splitext(os.path.basename(volume))[0]
        new_img_path = join(outdir, lm_basename + FILE_FORMAT)

        cmd = [
            'transformix',
            '-in', volume,
            '-tp', tform,
            '-out', outdir
        ]

        if threads:
            cmd.extend(['-threads', str(threads)])
        try:
            subprocess.check_output(cmd)
        except Exception as e:
            print('transformix failed inverting volume: {} Is transformix installed?. Error: {}'.format(volume, e))
            print(e)
            #logging.error('transformix failed with this command: {}\nerror message:'.format(cmd), exc_info=True)
            sys.exit()
        try:
            old_img = os.path.join(outdir, TRANSFORMIX_OUT)
            os.rename(old_img, new_img_path)
        except OSError:

            return old_img
        else:
            return new_img_path

