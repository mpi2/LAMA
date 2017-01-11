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
import subprocess as sub

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
    def __init__(self, wt_config_path, mut_proj_dir, wt_subset_file=None, mut_subset_file=None, in_dir=None, n1=True):
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

        if not is_r_installed():
            logging.error("Could not find an R installation. It is required for statistical analysis")
            sys.exit()

        # The root of the project dir for phenotype detection results
        self.mut_proj_dir = os.path.abspath(mut_proj_dir)

        # bool: do one against many analysis?
        self.n1 = n1

        self.wt_subset_file = wt_subset_file
        self.mut_subset_file = mut_subset_file

        if in_dir:
            self.in_dir = in_dir
        else:
            self.in_dir = join(self.mut_proj_dir, 'inputs')
            if not os.path.isdir(self.in_dir):
                print "\nCannot find input directory: {}\nEither place a file called inputs in project directory. Or specify input directory with -i option".format(self.in_dir)
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

    def write_config(self):
        """
        After getting the wildtype registration config and substituting mutant-specific info, write out a mutant
        registration config file
        """

        # Required parameters
        replacements = {'inputs': self.mut_config['inputs'],
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

        wt_intensity_abspath = self.wt_paths.get(INTENSITY_DIR)
        mut_intensity_abspath = self.mut_paths.get(INTENSITY_DIR)

        wt_intensity_dir = relpath(wt_intensity_abspath, stats_dir)
        mut_intensity_dir = relpath(mut_intensity_abspath, stats_dir)
        intensity_normalisation_roi = self.wt_config.get('normalisation_roi')

        # If there is no intensity directory, it means no normalization has occured. In this case, use the last
        # registration output directory

        if not os.path.exists(wt_intensity_abspath):
            wt_int = self.get_last_reg_stage(self.wt_paths)
            wt_intensity_dir = relpath(wt_int, stats_dir)
            mut_int = self.get_last_reg_stage(self.mut_paths)
            mut_intensity_dir = relpath(mut_int, stats_dir)

        wt_glcm_dir = relpath(join(self.wt_paths.get(GLCM_DIR)), stats_dir)
        mut_glcm_dir = relpath(join(self.mut_paths.get(GLCM_DIR)), stats_dir)

        wt_deformation_dir = relpath(join(self.wt_paths.get(DEFORMATION_DIR)), stats_dir)
        mut_deformation_dir = relpath(join(self.mut_paths.get(DEFORMATION_DIR)), stats_dir)

        # Jacobians and deformations can be generated at differnt scales eg 192-8, or 64-8
        # Just do stats for the first scale in the list. Can craft a stats.yaml by hand if multiple scales
        # should be analysed
        def_config_entry = self.mut_config.get('generate_deformation_fields')
        if def_config_entry:
            first_def_scale = def_config_entry.keys()[0]

        wt_jacobian_dir = relpath(join(self.wt_paths.get(JACOBIAN_DIR), first_def_scale), stats_dir)
        mut_jacobian_dir = relpath(join(self.mut_paths.get(JACOBIAN_DIR), first_def_scale), stats_dir)

        mut_inverted_labels = relpath(self.mut_paths.get('inverted_labels'), stats_dir)
        wt_inverted_labels = relpath(self.wt_paths.get('inverted_labels'), stats_dir)

        fixed_mask = relpath(join(self.wt_config_dir, self.wt_config['fixed_mask']), stats_dir)

        stats_meta_path = join(stats_dir, STATS_METADATA_PATH)

        skip_inv = self.mut_config.get('skip_transform_inversion')
        if skip_inv:
            logging.info('Skipping inversion of transforms')
            inverted_tform_config = None
        else:
            inverted_mut_tform_dir = self.mut_paths.get('inverted_transforms')
            inverted_tform_config = relpath(join(inverted_mut_tform_dir, 'invert.yaml'), stats_dir)

        voxel_size = self.mut_config.get('voxel_size')

        project_name = self.mut_config.get('project_name')
        if not project_name:
            project_name = '_'

        label_map_path = self.mut_config.get('label_map')
        label_map_path = relpath(join(label_map_path, stats_dir))
        organ_names = self.mut_config.get('organ_names')
        organ_names = relpath(join(organ_names, stats_dir))

        # Create a config file for the stats module to use
        stats_config_dict = {
            'project_name': project_name,
            'fixed_mask': fixed_mask,
            'n1': self.n1,
            'data': {
                'intensity':
                    {
                     'wt': wt_intensity_dir,
                     'mut': mut_intensity_dir,
                     'tests': list(stats_tests_to_perform),  # copy or we end up with a reference to the orignal in yaml
                     'normalisation_roi': intensity_normalisation_roi
                     },
                'jacobians':
                    {
                     'wt': wt_jacobian_dir,
                     'mut': mut_jacobian_dir,
                     'tests': list(stats_tests_to_perform)
                    },
                # 'glcm':
                #     {
                #      'wt': wt_glcm_dir,
                #      'mut': mut_glcm_dir,
                #      'tests': list(stats_tests_to_perform)
                #      },
                'organvolumes':
                    {
                     'wt': wt_inverted_labels,
                     'mut': mut_inverted_labels,

                     'tests': list(stats_tests_to_perform) # This is ignored at the moment but needs to be there
                     },
            },
            'invert_config_file': inverted_tform_config,
            'wt_groups': wt_groups_relpath,
            'mut_groups': mut_groups_relpath,
            'formulas': list(formulas),
            'voxel_size': voxel_size,
            'organ_names': organ_names,
            'label_map_path': label_map_path
        }

        # Add the jacobians and deformations based on how many scales were looked at

        if self.wt_subset_file:
            wt_subset_file = join(self.mut_config_dir, self.wt_subset_file)
            wt_subset_relpath = relpath(wt_subset_file, stats_dir)
            stats_config_dict['wt_subset_file'] = wt_subset_relpath

        if self.mut_subset_file:
            mut_subset_file = join(self.mut_config_dir, self.mut_subset_file)
            mut_subset_relpath = relpath(mut_subset_file, stats_dir)
            stats_config_dict['mut_subset_file'] = mut_subset_relpath

        common.mkdir_if_not_exists(stats_dir)

        # Look for groups file in the config dir and merge with the groups file for the mutants
        if self.wt_config.get('generate_deformation_fields'):
            for deformation_id in self.wt_config['generate_deformation_fields'].keys():
                deformation_id = str(deformation_id)  # yaml interperets numbers with underscores as ints
                wt_deformation_scale_dir = join(wt_deformation_dir, deformation_id)
                wt_jacobian_scale_dir = join(wt_jacobian_dir, deformation_id)
                mut_deformation_scale_dir = join(mut_deformation_dir, deformation_id)
                mut_jacobian_scale_dir = join(mut_jacobian_dir, deformation_id)

                jacobians_scale_config = {
                    'wt': wt_jacobian_scale_dir,
                    'mut': mut_jacobian_scale_dir,
                    'tests': list(stats_tests_to_perform)
                }
                # deformations_scale_config = {
                #     'wt': wt_deformation_scale_dir,
                #     'mut':mut_deformation_scale_dir,
                #     'tests': list(stats_tests_to_perform)
                # }
                #
                # stats_config_dict['data']['deformations_' + deformation_id] = deformations_scale_config
                stats_config_dict['data']['jacobians_' + deformation_id] = jacobians_scale_config


        with open(stats_meta_path, 'w') as fh:
            fh.write(yaml.dump(stats_config_dict, default_flow_style=False))
        return stats_meta_path

    def run_registration(self, config):
        l = lama.RegistraionPipeline(config, create_modified_config=False)

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
        mutant_config['inputs'] = mut_inputs_relpath

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


def is_r_installed():
    installed = True

    FNULL = open(os.devnull, 'w')
    try:
        sub.call(['Rscript'], stdout=FNULL, stderr=sub.STDOUT)
    except sub.CalledProcessError:
        installed = False
    except OSError:
        installed = False
        logging.warn('R or Rscript not installed. Will not be able to use linear model')
    return installed


if __name__ == '__main__':

    # Log all uncaught exceptions
    sys.excepthook = common.excepthook_overide

    import argparse

    parser = argparse.ArgumentParser("MRC Harwell registration pipeline")
    parser.add_argument('-c', '--config', dest='wt_config', help='Config file of the wildtype run (YAML format)', required=True)
    parser.add_argument('-i', '--input', dest='in_dir', help='directory containing input volumes', required=False)
    parser.add_argument('-p', '--proj-dir', dest='mut_proj_dir', help='directory to put results', required=True)
    parser.add_argument('-n1', '--specimen_n=1', dest='n1', help='Do one mutant against many wts analysis?', default=False)
    parser.add_argument('-wt_list', '--wildtpe_list', dest='wt_list', help='List of volume names that defines a subset of wt volumes to use', default=False)
    parser.add_argument('-mut_list', '--mutant_list', dest='mut_list', help='List of volume names that defines a subset of mut volumes to use', default=False)
    args, _ = parser.parse_known_args()
    PhenoDetect(args.wt_config, args.mut_proj_dir, args.wt_list, args.mut_list, args.in_dir)

    args = parser.parse_args()




