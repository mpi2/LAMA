#!/usr/bin/env python


"""pheno_detect.py

Phenotype detection component of the MRC Harwell registration pipeline.
Uses a config file that was used to perform registration on a series of wild type embryos. The location of the wiltype
data needed for analysis is determined from this config file.

Example
-------

    $ pheno_detect.py -c wt_registration_config_file_path -i dir_with_input_volumes -p project_directory_path

Notes
-----


"""

import os
from os.path import join, relpath
import copy
import logging
import shutil
import yaml
import sys
import lama
from stats.lama_stats import LamaStats
import common
from paths import RegPaths

VOLUME_CALCULATIONS_FILENAME = "organvolumes.csv"

MUTANT_CONFIG = 'mutant_config_modified.yaml'
"""str: Location to save the genrated config file for registering the mutants"""

INTENSITY_DIR = 'intensity'
"""str: directory to save the normalised registered images to"""

JACOBIAN_DIR = 'jacobians'
"""str: directory to save the determinant of the jacobians to"""

DEFORMATION_DIR = 'deformations'
"""str: directory to save the deformation fileds to"""

GLCM_DIR = 'glcms'

STATS_METADATA_HEADER = "This file can be run like: reg_stats.py -c stats.yaml"
STATS_METADATA_PATH = 'stats.yaml'

LOG_FILE = 'LAMA.log'

ORGAN_VOLS_OUT = 'organ_volumes.csv'


class PhenoDetect(object):
    """
    Phenotype detection

    """
    def __init__(self, wt_config_path, mut_proj_dir, in_dir=None, n1=True):
        """
        Parameters
        ----------
        wt_config_path: str
            path to a wildtype cofig file.
        mut_proj_dir: str
            path to root of project directory in which to save phenotype detection results
        in_dir: str
            path to directory containing mutant volumes to analsye
        n1: bool
            whether to perform optional one against many analysis
        """

        # The root of the project dir for phenotype detection results
        self.mut_proj_dir = os.path.abspath(mut_proj_dir)

        # bool: do one against many analysis?
        self.n1 = n1

        if in_dir:
            self.in_dir = in_dir
        else:
            self.in_dir = join(self.mut_proj_dir, 'inputs')
            if not os.path.isdir(self.in_dir):
                print "\nCannot find input directory: {}\nEither place a file called "
                "inputs in project directory. Or specify input directory with -i option".format(self.in_dir)
                sys.exit()

        (self.wt_config, self.wt_config_dir,
         self.mut_config, self.mut_config_dir) = self.get_config(wt_config_path, self.in_dir)

        self.wt_paths = RegPaths(self.wt_config_dir, self.wt_config)
        self.mut_paths = RegPaths(self.mut_config_dir, self.mut_config)

        self.out_dir = self.mut_paths.default_outdir
        self.mut_config['output_dir'] = self.out_dir

        self.mut_config_path = join(self.mut_proj_dir, MUTANT_CONFIG)
        self.write_config()

        if not self.mut_config.get('fixed_mask'):
            self.fixed_mask = None
            print 'WT fixed mask not present. Optimal results will not be obtained'
            # TODO: fix logging. Maybe put in root dir
            logging.warn('WT fixed mask not present. Optimal stats results might not be obtained')
        else:
            self.fixed_mask = self.mut_config['fixed_mask']

        self.run_registration(self.mut_config_path)
        # logfile = join(self.out_dir, LOG_FILE)
        # common.init_logging(logfile)

        logging.info('Stats analysis started')

        stats_metadata_path = self.write_stats_config()
        LamaStats(stats_metadata_path)

        # #Calculate organ volumes and do some stats on them
        # wt_organ_vols = join(self.wt_config_dir, self.wt_config['label_inversion_dir'],
        #                      VOLUME_CALCULATIONS_FILENAME)
        #
        # mut_organ_volumes = join(self.out_dir,
        #                          self.mut_config['label_inversion_dir'],
        #                          VOLUME_CALCULATIONS_FILENAME)
        # organ_vol_stats_out = join(self.out_dir, 'organ_volume_stats.csv')
        #
        # #organvolume_stats(wt_organ_vols, mut_organ_volumes, organ_vol_stats_out)

    def write_config(self):
        """
        After getting the wildtype registration config and substituting mutant-specific info, write out a mutant
        registration config file
        """

        # Required parameters
        replacements = {'inputvolumes_dir': self.mut_config['inputvolumes_dir'],
                        'fixed_volume': self.mut_config['fixed_volume'],
                        'fixed_mask': self.mut_config['fixed_mask'],
                        'pad_dims': self.mut_config['pad_dims']
                        }

        # optional parameters
        if self.mut_config.get('label_map_path'):
            replacements['label_map_path'] = self.mut_config['label_map_path']

        if self.mut_config.get('organ_names'):
            replacements['organ_names'] = self.mut_config['organ_names']

        if self.mut_config.get('isosurface_dir'):
            replacements['isosurface_dir'] = self.mut_config['isosurface_dir']

        lama.replace_config_lines(self.mut_config_path, replacements)

    def write_stats_config(self):
        """
        Writes a yaml config file for use by the reg_stats.py module. Provides paths to data and some options
        """
        stats_dir = self.mut_paths.get('stats')

        if self.wt_config.get('stats_tests'):
            stats_tests_to_perform = self.wt_config['stats_tests']
        else:
            stats_tests_to_perform = ['LM']
            logging.info("No stats test specified. Using default of Linear model")

        if self.wt_config.get('formulas'):
            formulas = self.wt_config['formulas']
        else:
            formulas = ['data ~ genotype'] # TODO. Just get the second column header name from groups.csv
            logging.info("No formulas specified. Using default: 'data ~ genotype")

        # Get the groups file paths
        wt_groups_name =  self.wt_config.get('groups_file')
        mut_groups_name = self.mut_config.get('groups_file')

        if all((wt_groups_name, mut_groups_name)):
            wt_groups = join(self.wt_config_dir, wt_groups_name)
            mut_groups = join(self.mut_config_dir, mut_groups_name)

        else:  # look for default groups file location
            wt_groups = join(self.wt_config_dir, 'groups.csv')
            mut_groups = join(self.mut_config_dir, 'groups.csv')

        # Check that the paths exists
        if all([os.path.exists(x) for x in [wt_groups, mut_groups]]):
            wt_groups_relpath = relpath(wt_groups, stats_dir)
            mut_groups_relpath = relpath(mut_groups, stats_dir)
        else:
            wt_groups_relpath = None
            mut_groups_relpath = None

        wt_intensity_dir = relpath(join(self.wt_paths.get(INTENSITY_DIR)), stats_dir)
        mut_intensity_dir = relpath(join(self.mut_paths.get(INTENSITY_DIR)), stats_dir)

        # If there is no intensity directory, it means no normalization has occured. In this case, use the last
        # registration output directory
        if not os.path.exists(wt_intensity_dir):
            wt_int = self.get_last_reg_stage(self.wt_paths)
            wt_intensity_dir = relpath(wt_int, stats_dir)
            mut_int = self.get_last_reg_stage(self.mut_paths)
            mut_intensity_dir = relpath(mut_int, stats_dir)

        wt_glcm_dir = relpath(join(self.wt_paths.get(GLCM_DIR)), stats_dir)
        mut_glcm_dir = relpath(join(self.mut_paths.get(GLCM_DIR)), stats_dir)

        wt_deformation_dir = relpath(join(self.wt_paths.get(DEFORMATION_DIR)), stats_dir)
        mut_deformation_dir = relpath(join(self.mut_paths.get(DEFORMATION_DIR)), stats_dir)

        wt_jacobian_dir = relpath(join(self.wt_paths.get(JACOBIAN_DIR)), stats_dir)
        mut_jacobian_dir = relpath(join(self.mut_paths.get(JACOBIAN_DIR)), stats_dir)

        mut_organ_vols_file = relpath(self.mut_paths.get(ORGAN_VOLS_OUT), stats_dir)
        wt_organ_vols_file = relpath(self.wt_paths.get(ORGAN_VOLS_OUT), stats_dir)

        fixed_mask = relpath(join(self.wt_config_dir, self.wt_config['fixed_mask']), stats_dir)

        stats_meta_path = join(stats_dir, STATS_METADATA_PATH)

        inv = self.mut_config.get('skip_transform_inversion')
        if inv:
            logging.info('Skipping inversion of transforms')
            inverted_tform_config = None
        else:
            inverted_mut_tform_dir = self.mut_paths.get('inverted_transforms')
            inverted_tform_config = relpath(join(inverted_mut_tform_dir, 'invert.yaml'), stats_dir)

        voxel_size = self.mut_config.get('voxel_size')

        project_name = self.mut_config.get('project_name')
        if not project_name:
            project_name = '_'

        # Save the path to the project log path in case we need to run stats seperately
        project_root = relpath(self.mut_config_dir, stats_dir)


        # Create a metadata file for the stats module to use
        stats_metadata = {
            'project_name': project_name,
            'project_root': project_root,
            'fixed_mask': fixed_mask,
            'n1': self.n1,
            'data': {
                'intensity':
                    {
                     'wt': wt_intensity_dir,
                     'mut': mut_intensity_dir,
                     'tests': list(stats_tests_to_perform)  # copy or we end up with a reference to the orignal in yaml
                     },
                # 'glcm':
                #     {
                #      'wt': wt_glcm_dir,
                #      'mut': mut_glcm_dir,
                #      'tests': list(stats_tests_to_perform)
                #      },
                'deformations':
                    {
                     'wt': wt_deformation_dir,
                     'mut': mut_deformation_dir,
                     'tests': list(stats_tests_to_perform)

                     },
                'jacobians':
                    {
                     'wt': wt_jacobian_dir,
                     'mut': mut_jacobian_dir,
                     'tests': list(stats_tests_to_perform)
                     },
                # 'organ_volumes':
                #     {
                #      'wt': wt_organ_vols_file,
                #      'mut': mut_organ_vols_file,
                #      'tests': list(stats_tests_to_perform)
                #      },
            },
            'invert_config_file': inverted_tform_config,
            'wt_groups': wt_groups_relpath,
            'mut_groups': mut_groups_relpath,
            'formulas': list(formulas),
            'voxel_size': voxel_size
        }

        common.mkdir_if_not_exists(stats_dir)

        # Look for groups file in the config dir and merge with the groups file for the mutants

        with open(stats_meta_path, 'w') as fh:
            fh.write(yaml.dump(stats_metadata, default_flow_style=False))
        return stats_meta_path

    def run_registration(self, config):
        reg = lama.RegistraionPipeline(config, create_modified_config=False)

    def get_config(self, wt_config_path, mut_in_dir):
        """
        Gets the config file that was used for the wildtype registration.
        Copies it and fills out the mutant-specific entries
        """
        wt_config_dir = os.path.abspath(os.path.dirname(wt_config_path))
        mut_config_path = join(self.mut_proj_dir, MUTANT_CONFIG)

        # Copy the wildtype config file to the mutant directory
        shutil.copyfile(wt_config_path, mut_config_path)

        mut_config_dir = os.path.abspath(os.path.dirname(mut_config_path))

        wt_config = yaml.load(open(wt_config_path, 'r'))
        mutant_config = copy.deepcopy(wt_config)

        mut_inputs_relpath = relpath(mut_in_dir, os.path.dirname(mut_config_path))
        mutant_config['inputvolumes_dir'] = mut_inputs_relpath

        fixed_vol_path = join(wt_config_dir, wt_config['fixed_volume'])
        fixed_vol_rel = relpath(fixed_vol_path, mut_config_dir)

        fixed_mask_path = join(wt_config_dir, wt_config['fixed_mask'])
        fixed_mask_rel = relpath(fixed_mask_path, mut_config_dir)

        if wt_config.get('label_map_path'):  # optional parameter
            fixed_labelmap = join(wt_config_dir, wt_config['label_map_path'])
            fixed_labelmap_rel = relpath(fixed_labelmap, mut_config_dir)
            mutant_config['label_map_path'] = fixed_labelmap_rel

        if wt_config.get('organ_names'):  # Optional parameter
            fixed_organ_names = join(wt_config_dir, wt_config['organ_names'])
            organ_names_rel = relpath(fixed_organ_names, mut_config_dir)
            mutant_config['organ_names'] = organ_names_rel

        if wt_config.get('isosurface_dir'):
            mesh_dir = join(wt_config_dir, wt_config['isosurface_dir'])
            mesh_dir_rel = relpath(mesh_dir, mut_config_dir)
            mutant_config['isosurface_dir'] = mesh_dir_rel

        mutant_config['fixed_volume'] = fixed_vol_rel
        mutant_config['fixed_mask'] = fixed_mask_rel
        mutant_config['pad_dims'] = wt_config['pad_dims']

        return wt_config, wt_config_dir, mutant_config, mut_config_dir

    def get_last_reg_stage(self, path_object):
        last_reg_id = self.wt_config['registration_stage_params'][-1]['stage_id']
        reg_path = join(path_object.get('root_reg_dir'), last_reg_id)

        return reg_path


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser("MRC Harwell registration pipeline")
    parser.add_argument('-c', '--config', dest='wt_config', help='Config file of the wildtype run (YAML format)', required=True)
    parser.add_argument('-i', '--input', dest='in_dir', help='directory containing input volumes', required=False)
    parser.add_argument('-p', '--proj-dir', dest='mut_proj_dir', help='directory to put results', required=True)
    parser.add_argument('-n1', '--specimen_n=1', dest='n1', help='Do one mutant against many wts analysis?', default=False)
    args, _ = parser.parse_known_args()
    PhenoDetect(args.wt_config, args.mut_proj_dir, args.in_dir)


    args = parser.parse_args()




