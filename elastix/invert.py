#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""invert.py

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

IGNORE_FOLDER = 'resolution_images'

import logging
import tempfile
import os
import subprocess
import sys
from collections import defaultdict
from multiprocessing import Pool
from os.path import join, splitext, abspath, basename, isfile, isdir
import shutil

import yaml
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common
from img_processing.pad import unpad_roi
from paths import RegPaths

ELX_TRANSFORM_PREFIX = 'TransformParameters.0.txt'
ELX_PARAM_PREFIX = 'elastix_params_'
ELX_INVERTED_POINTS_NAME = 'outputpoints.vtk'
FILE_FORMAT = '.nrrd'
LOG_FILE = 'inversion.log'
TRANSFORMIX_OUT = 'result.nrrd'
INVERSION_DIR_NAME = 'Inverted_transform_parameters'
LABEL_INVERTED_TRANFORM = 'labelInvertedTransform.txt'
IMAGE_INVERTED_TRANSFORM = 'ImageInvertedTransform.txt'
VOLUME_CALCULATIONS_FILENAME = "organvolumes.csv"

def setup_logging(outdir, logname):
    """
    If this module is being run directly from command line (ie. not from lama.py) setup logging to a new file

    Parameters
    ----------
    outdir: str
        directory to save log file in
    logname: str
        name of log file
    """

    if __name__ == '__main__':
        logpath = join(outdir, logname)
        logging.basicConfig(filename=logpath, level=logging.DEBUG,
                            format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p')


def batch_invert_transform_parameters(config_file, invert_config_file, outdir, threads=None, noclobber=False):
    """
    Create new elastix TransformParameter files that can then be used by transformix to invert labelmaps, stats etc

    Parameters
    ----------
    config_file: str
        path to original reg pipeline config file

    outdir: str
        Absolute path to output dir

    invert_config_file: str
        path to output inversion config to
    """
    common.test_installation('elastix')

    setup_logging(outdir, 'invert_transforms.log')

    with open(config_file, 'r') as yf:
        config = yaml.load(yf)
        config_dir = os.path.abspath(os.path.dirname(config_file))

    reg_dirs = get_reg_dirs(config, config_dir)

    # Get the image basename from the first stage registration folder (rigid?)
    first_stage = join(config_dir, reg_dirs[0])
    volume_names = [basename(x) for x in common.get_file_paths(first_stage, ignore_folder=IGNORE_FOLDER)]

    common.mkdir_if_not_exists(outdir)
    stages_to_invert = defaultdict(list)

    jobs = []
    if not threads:
        threads = 1
    else:
        threads = int(threads)

    for i, vol_name in enumerate(volume_names):

        label_replacements ={
            'FinalBSplineInterpolationOrder': '0',
            'FixedInternalImagePixelType': 'short',
            'MovingInternalImagePixelType':  'short',
            'ResultImagePixelType': 'unsigned char',
            'WriteTransformParametersEachResolution': 'false',
            'WriteResultImageAfterEachResolution': 'false'
        }

        image_replacements = {
            'FinalBSplineInterpolationOrder': '3',
            'FixedInternalImagePixelType': 'float',
            'MovingInternalImagePixelType':  'float',
            'ResultImagePixelType': 'float',
            'WriteTransformParametersEachResolution': 'false',
            'WriteResultImageAfterEachResolution': 'false'

        }

        vol_id, vol_ext = splitext(vol_name)

        for r in reg_dirs:
            reg_dir = join(config_dir, r)

            stage_out_dir = join(outdir, basename(reg_dir))

            moving_dir = join(config_dir, reg_dir, vol_id)
            invert_param_dir = join(stage_out_dir, vol_id)

            if not os.path.isdir(moving_dir):
                logging.warning('cannot find {}'.format(moving_dir))
                continue
            stage_vol_files = os.listdir(moving_dir)  # All the files from the registration dir
            stage_files = os.listdir(reg_dir)         # The registration stage parent directory
            parameter_file = next(join(reg_dir, i) for i in stage_files if i.startswith(ELX_PARAM_PREFIX))
            transform_file = next(join(moving_dir, i) for i in stage_vol_files if i.startswith(ELX_TRANSFORM_PREFIX))

            if not isfile(parameter_file):
                logging.error('elastix transform parameter file missing: {}'.fomrat(transform_file))
                continue
            if not isfile(parameter_file):
                logging.error('elastix registration paramter file missing: {}'.format(parameter_file))
                continue

            common.mkdir_if_not_exists(stage_out_dir)

            rel_inversion_path = os.path.basename(r)
            if rel_inversion_path not in stages_to_invert['inversion_order']:
                stages_to_invert['inversion_order'].insert(0, rel_inversion_path)

            if noclobber:
                if os.path.isdir(invert_param_dir):
                    continue

            common.mkdir_force(invert_param_dir)  # Overwrite any inversion file that exist for a single specimen
            reg_metadata = yaml.load(open(join(moving_dir, common.INDV_REG_METADATA)))
            fixed_volume = join(moving_dir, reg_metadata['fixed_vol'])  # The original fixed volume used in the registration

            # Invert the Transform paramteres with options for normal image inversion

            job = {
                'invert_param_dir': invert_param_dir,
                'parameter_file': abspath(parameter_file),
                'transform_file': transform_file,
                'fixed_volume': fixed_volume,
                'param_file_output_name': 'imageParam.txt',
                'image_replacements': image_replacements,
                'label_replacements': label_replacements,
                'image_transform_file': IMAGE_INVERTED_TRANSFORM,
                'label_transform_file': LABEL_INVERTED_TRANFORM
            }

            jobs.append(job)

    # with open('/home/neil/work/jobs.json', 'w') as fh:
    #     json.dump(jobs, fh, sort_keys=True, indent=4, separators=(',', ': '))
    # return
    logging.info('inverting with {} threads: '.format(threads))
    pool = Pool(threads)
    try:
        pool.imap_unordered(_invert_transform_parameters, jobs)

    except KeyboardInterrupt:
        print('terminating inversion')
        pool.terminate()
        pool.join()

    reg_dir = os.path.relpath(reg_dir, outdir)
    stages_to_invert['registration_directory'] = reg_dir
    # Create a yaml config file so that inversions can be run seperatley
    with open(invert_config_file, 'w') as yf:
        yf.write(yaml.dump(dict(stages_to_invert), default_flow_style=False))


def _invert_transform_parameters(args):
    """
    Generate a single inverted elastix transform parameter file. This can then be used to invert labels, masks etc.
    If any of the step faile, return as subsequent steps will also fail. The logging of failrures is handled
    within each function
    """
    label_param = abspath(join(args['invert_param_dir'], args['param_file_output_name']))
    # Modify the elastix registration parameter file with the image parameters
    _modify_param_file(abspath(args['parameter_file']), label_param, args['image_replacements'])

     # Make the inverted TransformParameters file
    if not _invert_tform(args['fixed_volume'], abspath(args['transform_file']), label_param, args['invert_param_dir']):
        return

    # Get the resulting TransformParameters file and add initial transforms if needed
    image_inverted_tform = abspath(join(args['invert_param_dir'], 'TransformParameters.0.txt'))
    image_transform_param_path = abspath(join(args['invert_param_dir'], args['image_transform_file']))
    if not _modify_tform_file(image_inverted_tform, image_transform_param_path):
        return

    # Now create a TransformParameter file that is suitable for inverting label maps and masks (interpolation 0)
    label_transform_param_path = join(abspath(join(args['invert_param_dir'], args['label_transform_file'])))
    # replace the parameter in the image file with label-specific parameters and save in new file. No need to
    # generate one from scratch
    if not _modify_param_file(image_transform_param_path, label_transform_param_path, args['label_replacements']):
        return

    _modify_tform_file(label_transform_param_path)



def get_reg_dirs(config, config_dir):
    """

    """
    paths = RegPaths(config_dir, config)
    reg_stages = []
    root_reg_dir = paths.get('root_reg_dir')
    for i, reg_stage in enumerate(config['registration_stage_params']):
        stage_id = reg_stage['stage_id']
        stage_dir = join(root_reg_dir, stage_id)
        reg_stages.append(stage_dir)
    return reg_stages


class Invert(object):
    def __init__(self, config_path, invertables, outdir, threads=None):
        """
        Inverts a of volumes. A yaml config file specifies the order of inverted transform parameters
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

        :return:
        """

        with open(config_path, 'r') as yf:
            self.config = yaml.load(yf)

        self.invertables = invertables
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

        inverting_names = os.listdir(self.inverted_tform_stage_dirs[0])

        for i, vol_name in enumerate(inverting_names):
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

                invertable = self._invert(invertable, transform_file, invert_vol_out_dir, self.threads)
        self.last_invert_dir = invert_stage_out

    def _invert(self):
        raise NotImplementedError

    @staticmethod
    def _rename_output(output, folder_name):
        """
        Parameters
        ----------
        output: str
            path to the file to change name
        folder_name: str
            path to output folder. The last bit of this path contains name to convert to
        """
        path, base = os.path.split(os.path.normpath(folder_name))
        new_output_name = os.path.join(folder_name, '{}.nrrd'.format(base)) # as we ar enot prepending with Seq we can just do a move now?
        try:
            os.rename(output, new_output_name)
        except OSError:
            print 'could not rename {}'.format(output)
            return output
        else:
            return new_output_name


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
        common.test_installation('transformix')
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
            sys.exit()
            #logging.error('transformix failed with this command: {}\nerror message:'.format(cmd), exc_info=True)

        old_img = os.path.join(outdir, TRANSFORMIX_OUT)
        new_output_name = self._rename_output(old_img, outdir)

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

    def __init__(self,  config_path, invertables, outdir, threads=None):
        super(InvertMeshes, self).__init__(config_path, invertables, outdir, threads)
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
            print 'transformix failed inverting mesh: {}'.format(mesh)
            logging.error('transformix failed with this command: {}\nerror message:'.format(cmd), exc_info=True)
            print e
            sys.exit(1)
        try:
            # rename the inverted points form this stage
            old_vtk = os.path.join(outdir, ELX_INVERTED_POINTS_NAME)
            os.rename(old_vtk, new_vtk_path)
        except OSError:

            raise
        else:
            return new_vtk_path


class InvertRoi(InvertLabelMap):
    def __init__(self, config_path, invertables, outdir,  vol_info, voxel_size, threads=None):
        super(InvertRoi, self).__init__(config_path, invertables, outdir, threads)
        self.invert_transform_name = LABEL_INVERTED_TRANFORM
        self.vol_info = vol_info
        self.voxel_size = voxel_size

    def run(self):
        super(InvertRoi, self).run()
        # At this point we have a bunch of rois inverted onto the padded inputs
        # We need to adjust the rois to account for the padding
        out = join(self.out_dir, 'Extracted_roi')
        unpad_roi(self.vol_info, self.last_invert_dir, self.voxel_size, out)


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
            print 'transformix failed inverting volume: {} Is transformix installed?. Error: {}'.format(volume, e)
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


def _modify_param_file(elx_param_file, newfile_name, replacements):
    """
    Modifies the elastix input parameter file that was used in the original transformation.
    Adds DisplacementMagnitudePenalty (which is needed for inverting)
    Turns off writing the image results at the end as we only need an inveterted output file.
    Also changes interpolation order in the case of inverting labels

    Parameters
    ----------
    elx_param_file: str
        path to elastix input parameter file
    newfile_name: str
        path to save modified parameter file to

    """

    try:
        with open(elx_param_file) as old, open(newfile_name, "w") as new:

            for line in old:
                if line.startswith("(Metric "):
                    line = '(Metric "DisplacementMagnitudePenalty")\n'
                if line.startswith('(WriteResultImage '):
                    line = '(WriteResultImage "false")\n'
                if line.startswith('WriteResultImageAfterEachResolution '):
                   continue
                try:
                    param_name = line.split()[0][1:]
                except IndexError:
                    continue  # comment?

                if param_name in replacements:
                    value = replacements[param_name]
                    try:
                        int(value)
                    except ValueError:
                        # Not an int, neeed quotes
                        line = '({} "{}")\n'.format(param_name, value)
                    else:
                        # An int, no quotes
                        line = '({} {})\n'.format(param_name, value)
                new.write(line)
    except IOError as e:
        logging.error("Error modifying the elastix parameter file: {}".format(e))
        return False
    return True


def _invert_tform(fixed, tform_file, param, outdir):
    """
    Invert the transform and get a new transform file
    """
    common.test_installation('elastix')

    cmd = ['elastix',
           '-t0', tform_file,
           '-p', param,
           '-f', fixed,
           '-m', fixed,
           '-out', outdir,
           '-threads', '1'
           ]

    # Just use one thread within elastix as LAMA is dealing with the multithreading
    try:
        subprocess.check_output(cmd)
    except (Exception, subprocess.CalledProcessError) as e:
        logging.exception('Inverting transform file failed. cmd: {}\n{}:'.format(cmd, str(e)))
          # We need to stop here as it can affect later stages. TODO: ned to make other threads halt too
        return False
    return True


def _modify_tform_file(elx_tform_file, newfile_name=None):
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

    if not newfile_name:  # Write to temporary file before overwriting
        new_file = tempfile.NamedTemporaryFile().name
    else:
        new_file = newfile_name

    new_tform_param_fh = open(new_file, "w+")

    try:
        tform_param_fh = open(elx_tform_file, "r")
        for line in tform_param_fh:
            if line.startswith('(InitialTransformParametersFileName'):
                line = '(InitialTransformParametersFileName "NoInitialTransform")\n'
            new_tform_param_fh.write(line)
        new_tform_param_fh.close()
        tform_param_fh.close()

    except IOError:
        logging.warning("Can't find tform file {}".format(elx_tform_file))
        return False

    if not newfile_name:
        shutil.move(new_file, elx_tform_file)

    return True


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


if __name__ == '__main__':

    # Log all uncaught exceptions
    sys.excepthook = common.excepthook_overide

    def print_args_error():
        msg = ("\nOptions are\n\n"
               "reg - make inverse transform parameter files for elastix\n"
               "labels - invert a label image using previously-generated inverse transform parameter files\n"
               "vol - invert a grey scale image using previously-generated inverse transform parameter files\n"
               "meshes - invert a itk mesh using previously-generated inverse transform parameter files\n"
               "roi - trnsform roi coordinates using previously-generated inverse transform parameter files\n\n"
               "Examples:\ninvert.py reg -c lama_config.yaml -o output/inverted_transforms, -t 8\n"
               "invert.py lables -c inverted_transforms/invert.yaml -o output/inverted_lables -i label_to_invert -t 8")
        sys.exit(msg)

    if len(sys.argv) < 2:
        print_args_error()

    import argparse

    # this_script_path = os.path.dirname(os.path.realpath(__file__))
    # default_invert_parameter = join(this_script_path, )

    if sys.argv[1] == 'labels':
        parser = argparse.ArgumentParser("invert elastix registrations and calculate organ volumes")
        parser.add_argument('-c', '--config', dest='config', help='yaml config file. Usually root/output/inverted_transforms/invert.yaml', required=True)
        parser.add_argument('-i', '--invertable', dest='invertable', help='label volume to invert', required=True)
        parser.add_argument('-o', '--outdir', dest='outdir', help='output dir. Usually root/output/inverted_labels', required=True)
        parser.add_argument('-t', '--threads', dest='threads', type=str, help='number of threads to use', required=False)

        args, _ = parser.parse_known_args()
        inv = InvertLabelMap(args.config, args.invertable, args.outdir, threads=args.threads)
        inv.run()

    elif sys.argv[1] == 'reg':
        parser = argparse.ArgumentParser("invert elastix registrations to create elastix inversion parameter files")
        parser.add_argument('-c', '--config',  dest='config', help='Main LAMA config file with list of registration dirs', required=True)
        parser.add_argument('-o', '--out',  dest='outdir', help='where to put the output', required=True)
        parser.add_argument('-t', '--threads', dest='threads', type=str, help='number of threads to use', required=False)
        parser.add_argument('-noclobber', '--noclobber', dest='noclobber', default=False, action='store_true')
        args, _ = parser.parse_known_args()
        config_out = join(args.outdir, 'invert.yaml')
        batch_invert_transform_parameters(args.config, config_out, args.outdir, args.threads, noclobber=args.noclobber)

    elif sys.argv[1] == 'vol':
        parser = argparse.ArgumentParser("invert elastix registrations and calculate organ volumes")
        parser.add_argument('-c', '--config', dest='config', help='yaml config file. Usually root/output/inverted_transforms/invert.yaml', required=True)
        parser.add_argument('-i', '--invertable', dest='invertable', help='volume to invert', required=True)
        parser.add_argument('-o', '--outdir', dest='outdir', help='output dir. Usually root/output/inverted_labels', required=True)
        parser.add_argument('-p', '--prefix', dest='prefix', help='A prefix added to the invertable, that is not present on the invert transform files', default=False)

        parser.add_argument('-t', '--threads', dest='threads', type=str, help='number of threads to use', required=False)
        args, _ = parser.parse_known_args()
        inv = InvertSingleVol(args.config, args.invertable, args.outdir)
        inv.run(args.prefix)

    elif sys.argv[1] == 'meshes':
        parser = argparse.ArgumentParser("invert meshes")
        parser.add_argument('-c', '--config', dest='config', help='yaml config file', required=True)
        parser.add_argument('-m', '--meshes', dest='mesh', help='mesh dir/mesh file', required=True)
        parser.add_argument('-o', '--outdir', dest='outdir', help='output dir', required=True)
        parser.add_argument('-t', '--threads', dest='threads', type=str, help='number of threads to use', required=False)

        args, _ = parser.parse_known_args()
        if os.path.isdir(args.mesh):
            for path in common.get_file_paths(args.mesh):
                inv = InvertMeshes(args.config, path, args.outdir)
                inv.run()
        else:
            inv = InvertMeshes(args.config, args.mesh, args.outdir)
            inv.run()

    elif sys.argv[1] == 'roi':
        parser = argparse.ArgumentParser("invert roi")
        parser.add_argument('-c', '--config', dest='config', help='yaml config file', required=True)
        parser.add_argument('-s', '--starts', dest='starts', help='roi starts (xyz)', required=True)
        parser.add_argument('-e', '--ends', dest='ends', help='roi ends (xyz)', required=True, nargs=3, type=int)
        parser.add_argument('-o', '--outdir', dest='outdir', help='output dir', required=True, nargs=3, type=int)
        parser.add_argument('-t', '--threads', dest='threads', type=str, help='number of threads to use', required=False)
        parser.add_argument('-i', '--info', dest='info', type=str, help='info on padding and full res locations, yaml',
                            required=False)
        parser.add_argument('-v', '--voxel_size', dest='voxel_size', type=str, help='Voxel size of scaled images (um)',
                            required=False)
        args, _ = parser.parse_known_args()
        inv = InvertRoi(args.config, args.label, args.outdir, args.info, args.voxel_size, args.threads)
        inv.run()

    else:
        print_args_error()

