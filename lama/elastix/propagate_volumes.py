#!/usr/bin/env python3

"""propagate_volumes.py

This module propagtes labels using transforms made from either 'invert_transforms.py' or 'reverse_registration.py'.

Note the word 'invert' used throughout. Read as propogate. It's this way as initially the only way to propagate
labels was to use the transform inversion method. But now we can use the reverse registration method i.e register
population average to specimen image thhen priopagate atlas using the given transform.

Example
-------

InvertLabelMap(invert_config, label_map_path, labels_inverion_dir, threads=32).run()

example config file:
    label_propagation_order:
    - rigid
    - affine
    - deformable_to_8
    - deformable_to_128

Note: The order will be reversed if doing label propagation using the inverted transform method

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
from typing import List
import os
import subprocess
from os.path import join
import shutil

from logzero import logger as logging
import yaml

from lama import common
from lama.common import cfg_load
from lama.elastix import (PROPAGATE_LABEL_TRANFORM, PROPAGATE_IMAGE_TRANSFORM, ELX_PARAM_PREFIX, TRANSFORMIX_OUT,
                          ELX_TRANSFORM_NAME, ELX_INVERTED_POINTS_NAME, PROPAGATE_CONFIG)


class Propagate(object):
    def __init__(self, config_path: Path, invertable, outdir, threads=None, noclobber=False):
        """
        Inverts a series of volumes. A yaml config file specifies the order of inverted transform parameters
        to use. This config file should be in the root of the directory containing these inverted tform dirs.

        Also need to input a directory containing volumes/label maps etc to invert. These need to be in directories
        named with the same name as the corresponding inverted tform file directories

        Parameters
        ----------
        config_path
            path to yaml config containing the oder of the inverted directories to use. The directories containing
            propagation tfransfrom files should be in the same diretory
        threads: str/ None
            number of threas to use. If None, use all available threads
        invertable: str
            path to object to invert (raw image, mask, label map etc)
        outdir
            where to store inverted volumes
        invertable: str
            dir or path. If dir, invert all objects within the subdirectories.
                If path to object (eg. labelmap) invert that instead
        noclobber: bool
            if True do not overwrite already inverted labels

        """

        self.noclobber = noclobber

        common.test_installation('transformix')

        self.config = cfg_load(config_path)

        self.invertables = invertable
        self.config_dir = config_path.parent  # The dir containing the inverted elx param files

        self.threads = threads
        self.out_dir = outdir
        common.mkdir_if_not_exists(self.out_dir)

        self.elx_param_prefix = ELX_PARAM_PREFIX
        self.PROPAGATION_TFORM_NAME = None  # Set in subclasses
        self.last_invert_dir = None # I thik this is used as a way to find volumes to do organ vol calculation on

    def run(self):
        """

        """
        # temp_tform_dir = Path('/home/neil/Desktop/t/test_propagarte/temp')
        done_file = self.out_dir / 'propagation.done'

        common.touch(done_file)

        vol_ids = os.listdir(self._transform_dirs()[0])

        logging.info('propagating volumes')

        for i, id_ in enumerate(vol_ids):

            prop_out_dir: Path = self.out_dir / id_  # create a folder with vol_id in case we have multiple vols to do
            prop_out_dir.mkdir(exist_ok=True)

            tform_root = self.config_dir
            init_tform = chain_tforms(tform_root, prop_out_dir, self.PROPAGATION_TFORM_NAME, self.config)

            propagated = self._propagate(self.invertables, init_tform, prop_out_dir, self.threads)

            if not propagated: # If inversion failed or there is nocobber, will get None
                continue

        self.last_invert_dir = self.out_dir

    def _propagate(self):
        raise NotImplementedError

    def _transform_dirs(self) -> List[Path]:
        """
        When implemented, this method returns a list of directories each containing transform parameter file for a
        label propagation stage.
        """
        raise NotImplementedError


class PropagateLabelMap(Propagate):

    def __init__(self, *args, **kwargs):
        super(PropagateLabelMap, self).__init__(*args, **kwargs)
        self.PROPAGATION_TFORM_NAME = PROPAGATE_LABEL_TRANFORM

    def _transform_dirs(self):
        dirs = []

        for dir_name in self.config['label_propagation_order']:
            dir_path = self.config_dir / dir_name
            dirs.append(dir_path)
        return dirs

    def _propagate(self, labelmap, tform, outdir: Path, threads=None):
        """
        Using the iverted elastix transform paramter file, invert a volume with transformix

        Parameters
        ----------
        vol: str
            path to volume to invert
        tform: str
            path to save inverted specimen volume
        outdir:
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

        old_img = outdir / TRANSFORMIX_OUT       # where the transformix-inverted labelmap will be

        id_ = outdir.name
        new_output_name = outdir / f'{id_}.nrrd'  # Renamed transformix-inverted labelmap

        # if self.noclobber and isfile(new_output_name):
        # Maybe need to do two types of noclobber
        # 1: where if the folder does not exist, do not do it
        # 2: where the folder exists but the final output file does not exist
        #     return None

        cmd = [
            'transformix',
            '-in', str(labelmap),
            '-tp', str(tform),
            '-out', str(outdir)
        ]

        if threads:
            cmd.extend(['-threads', str(threads)])
        try:
            subprocess.check_output(cmd)
        except Exception as e:
            logging.exception('{}\ntransformix failed propagating labelmap: {}'.format(e, labelmap))
            raise

        try:
            shutil.move(old_img, new_output_name)
        except IOError as e:
            print()
            'could not rename {}'.format(old_img)
            return old_img
        else:
            return new_output_name


class PropagateHeatmap(PropagateLabelMap):
    """
    This class behaves the same as InvertLabelap but uses a different transform parameter file
    """
    def __init__(self, *args, **kwargs):
        super(PropagateHeatmap, self).__init__(*args, **kwargs)
        self.invert_transform_name = PROPAGATE_IMAGE_TRANSFORM


class PropagateMeshes(Propagate):

    def __init__(self, config_path, invertable, outdir, threads=None):
        super(PropagateMeshes, self).__init__(config_path, invertable, outdir, threads)
        self.invert_transform_name = ELX_TRANSFORM_NAME

    def _transform_dirs(self):
        """
        For mesh inverion we need to use the orginal registration tranform files not the inverted ones

        Returns
        -------

        """
        dirs = []
        reg_dir = self.config.get('registration_directory')
        for dir_name in self.config['label_propagation_order']:
            dir_path = join(self.config_dir, reg_dir, dir_name)
            dirs.append(dir_path)
        return dirs

    def _propagate(self, mesh, tform, outdir, threads=None):
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
            cmd.extend(['-threads', str(threads)])
        try:
            subprocess.check_output(cmd)
        except Exception as e:
            logging.exception('transformix failed with this command: {}\nerror message:'.format(cmd), exc_info=True)
            raise
        try:
            # rename the inverted points form this stage
            old_vtk = os.path.join(outdir, ELX_INVERTED_POINTS_NAME)
            os.rename(old_vtk, new_vtk_path)
        except OSError:

            raise
        else:
            return new_vtk_path


class PropagateSingleVol(Propagate):
    """
    Invert volumes using the elastix inverted transform parameters.
    This class is used for inverting statitistics overlays
    """
    def __init__(self, *args, **kwargs):
        super(PropagateSingleVol, self).__init__(*args, **kwargs)
        self.invert_transform_name = PROPAGATE_IMAGE_TRANSFORM

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

            transform_file = join(inv_tform_dir, PROPAGATE_IMAGE_TRANSFORM)
            invert_vol_out_dir = join(invert_stage_out, volname)
            common.mkdir_if_not_exists(invert_vol_out_dir)

            invertable = self._propagate(invertable, transform_file, invert_vol_out_dir, self.threads)

    def _propagate(self, volume, tform, outdir, threads=None):
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
            logging.exception('transformix failed propagating volume: {} Is transformix installed?. Error: {}'.format(volume, e))
            raise
        try:
            old_img = os.path.join(outdir, TRANSFORMIX_OUT)
            os.rename(old_img, new_img_path)
        except OSError:

            return old_img
        else:
            return new_img_path


def chain_tforms(root_dir: Path, new_tform_dir, tform_name, config):

    label_replacements = {
        'FinalBSplineInterpolationOrder': '0',
        'FixedInternalImagePixelType': 'short',
        'MovingInternalImagePixelType': 'short',
        'ResultImagePixelType': '"unsigned char"',
        'WriteTransformParametersEachResolution': 'false',
        'WriteResultImageAfterEachResolution': 'false'
    }

    stages = config['label_propagation_order']

    # stages = stages[::-1]
    for i, stage in enumerate(stages):
        stage_dir = root_dir / stage
        tform_file = next(stage_dir.glob(f'**/{tform_name}'))
        new_tform_file = new_tform_dir /  f'{stage}_{tform_file.name}'
        shutil.copyfile(tform_file, new_tform_file)

        if i + 1 < len(stages):
            init_tform = new_tform_dir /  f'{stages[i+1]}_{tform_file.name}'
            # init_tform = f'{stages[i+1]}.txt'

        else:
            init_tform = None

        new_tform_path = new_tform_dir / f'{stage}_{tform_file.name}'

        if i == 0:
            file_for_transformix = new_tform_path

        with open(tform_file, 'r') as fh, open(new_tform_path, 'w') as nfh:

            for line in fh:
                #
                if line.startswith('(InitialTransformParametersFileName') and init_tform:
                    line = f'(InitialTransformParametersFileName "{str(init_tform)}")\n'
                else:
                    for param in label_replacements:
                        if line.startswith(f'({param}'):
                            line = f'({param} {label_replacements[param]})\n'

                nfh.write(line)

    return file_for_transformix