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
from os.path import join, relpath, abspath, exists
import copy
import logging
import shutil
import yaml
import sys
import lama
from stats.lama_stats import run as run_lama_stats
import common
from paths import RegPaths
import subprocess as sub
from lib import addict

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

DEFAULT_INPUT_DIR = 'inputs'

LOG_FILE = 'LAMA.log'

ORGAN_VOLS_OUT = 'organ_volumes.csv'


class PhenoDetect(object):
    """
    Phenotype detection

    """
    def __init__(self, wt_config_path, mut_proj_dir, wt_list_file=None, use_auto_staging=True, litter_baselines=None,
                 in_dir=None, n1=True):
        """
        Parameters
        ----------
        wt_config_path: str
            path to a wildtype cofig file.
        mut_proj_dir: str
            path to root of project directory in which to save phenotype detection results
        wt_list_file: str
            path to file with a list of volume basenames (not the full paths)of WTs to use in the analysis. Using this 
            option with prevent any automatic staging
        in_dir: str
            path to directory containing mutant volumes to analsye. deafualts to mut_proj_dir/inputs
        use_auto_staging: bool
            if auto staging information is available in baseline run, use this to select wild type test set
        n1: bool
            whether to perform optional one against many analysis
        litter_baselines: str
            path to folder containing litternate controls
        """

        if not common.is_r_installed():
            logging.error("Could not find an R installation. It is required for statistical analysis")
            sys.exit()

        # The root of the project dir for phenotype detection results
        self.mut_proj_dir = os.path.abspath(mut_proj_dir)

        # bool: do one against many analysis?
        self.n1 = n1

        self.wt_list_file = wt_list_file

        self.use_auto_staging = use_auto_staging

        self.litter_baselines_file = litter_baselines

        if in_dir:
            self.in_dir = in_dir
        else:
            self.in_dir = join(self.mut_proj_dir, DEFAULT_INPUT_DIR)
            if not os.path.isdir(self.in_dir):
                print "\nCannot find input directory: {}\nEither place a file called inputs in project directory. Or specify input directory with -i option".format(self.in_dir)
                sys.exit()

        (self.wt_config, self.wt_config_dir,
         self.mut_config, self.mut_config_dir) = self.get_config(wt_config_path, self.in_dir)

        self.wt_path_maker = RegPaths(self.wt_config_dir, self.wt_config)

        self.mut_path_maker = RegPaths(self.mut_config_dir, self.mut_config)

        self.out_dir = self.mut_path_maker.default_outdir
        self.mut_config['output_dir'] = self.out_dir

        self.wt_out_dir = self.wt_path_maker.default_outdir

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

        logging.info('Stats analysis started')

        stats_config_path = self.write_stats_config()
        run_lama_stats(stats_config_path)

    def write_config(self):
        """
        After getting the wildtype registration config and substituting mutant-specific info, write out a mutant
        registration config file. This will be used to tun lama on the mutants and will enable rerunning of the 
        process and act as a log
        """

        # Required parameters
        replacements = {'inputs': self.mut_config['inputs'],
                        'fixed_volume': self.mut_config['fixed_volume'],
                        'fixed_mask': self.mut_config['fixed_mask'],
                        'pad_dims': self.mut_config['pad_dims']
                        }

        # optional parameters
        if self.mut_config.get('label_map'):
            replacements['label_map'] = self.mut_config['label_map']

        if self.mut_config.get('organ_names'):
            replacements['organ_names'] = self.mut_config['organ_names']

        if self.mut_config.get('isosurface_dir'):
            replacements['isosurface_dir'] = self.mut_config['isosurface_dir']

        if self.mut_config.get('staging_volume'):
            replacements['staging_volume'] = self.mut_config['staging_volume']

        lama.replace_config_lines(self.mut_config_path, replacements)

    def write_stats_config(self):
        """
        Writes a yaml config file for use by the reg_stats.py module. Provides paths to data and some options
        bu_lama_stats.py will be run automatically asfter registration. But this config file will allow stats to be
        rerun at any time and to change parameters if required
        """
        stats_dir = self.mut_path_maker.get('stats')

        # stats_tests_to_perform = self.wt_config.get('stats_tests')  # Defaults to liner model (LM) in bu_lama_stats.py

        formulas = self.wt_config.get('formulas')  # Defaults to 'data~geneotype+crl' (crown-rump length) in lama_stats

        # Get the groups file paths for if we need to specifiy group membership in the linear model
        wt_groups_name = self.wt_config.get('groups_file')
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

        mut_inverted_labels = relpath(self.mut_path_maker.get('inverted_labels'), stats_dir)
        wt_inverted_labels = relpath(self.wt_path_maker.get('inverted_labels'), stats_dir)

        fixed_mask = relpath(join(self.wt_config_dir, self.wt_config['fixed_mask']), stats_dir)

        stats_meta_path = join(stats_dir, STATS_METADATA_PATH)

        skip_inv = self.mut_config.get('skip_transform_inversion')
        if skip_inv:
            logging.info('Skipping inversion of transforms')
            inverted_tform_config = None
        else:
            inverted_mut_tform_dir = self.mut_path_maker.get('inverted_transforms')
            inverted_tform_config = relpath(join(inverted_mut_tform_dir, 'invert.yaml'), stats_dir)

        voxel_size = self.mut_config.get('voxel_size')

        project_name = self.mut_config.get('project_name')

        label_map_path = self.mut_config.get('label_map')
        if label_map_path:
            label_map_path = relpath(join(self.wt_config_dir, self.wt_config['label_map']), stats_dir)

        organ_names = self.mut_config.get('label_names')
        if organ_names:
            organ_names = relpath(join(self.wt_config_dir, self.wt_config['label_names']), stats_dir)

        common.mkdir_if_not_exists(stats_dir)

        # Add the jacobians and deformations based on how many scales were looked at
        stats_config_dict = addict.Dict()
        self.add_intensity_stats_config(stats_config_dict, stats_dir)
        self.add_deformations_stats_config(stats_config_dict, stats_dir)

        if project_name:
            stats_config_dict['project_name'] = project_name

        if self.n1:
            stats_config_dict['n1'] = self.n1
        if formulas:
            stats_config_dict['formulas'] = formulas

        if label_map_path:
            stats_config_dict['label_map'] = label_map_path

        if wt_groups_relpath:
            stats_config_dict['wt_groups'] = wt_groups_relpath
        if mut_groups_relpath:
            stats_config_dict['mut_groups'] = mut_groups_relpath

        if inverted_tform_config: # If the run was inverted, add the invert config path so the stast can be inverted
            stats_config_dict['invert_config_file'] = inverted_tform_config
        if organ_names:
            stats_config_dict['organ_names'] = organ_names

        if voxel_size:
            stats_config_dict['voxel_size'] = voxel_size

        stats_config_dict['fixed_mask'] = fixed_mask

        if self.litter_baselines_file:
            stats_config_dict['littermate_controls'] = relpath(self.litter_baselines_file, stats_dir)

        if self.use_auto_staging:
            # get the paths to the mutant and wildtype staging files that were generated during
            wt_staging_file = join(self.wt_out_dir, common.STAGING_INFO_FILENAME)
            mut_staging_file = join(self.out_dir, common.STAGING_INFO_FILENAME)

            if not exists(wt_staging_file):
                logging.warn('cannot find wild type staging file {}'.format(wt_staging_file))
            if not exists(mut_staging_file):
                logging.warn('cannot find mutant type staging file {}'.format(mut_staging_file))
            if all([exists(wt_staging_file), exists(mut_staging_file)]):
                #  We have both staging files, add them to the stas config
                stats_config_dict['wt_staging_file'] = relpath(wt_staging_file, stats_dir)
                stats_config_dict['mut_staging_file'] = relpath(mut_staging_file, stats_dir)

        # Create intensity analysis section
        self.add_intensity_stats_config(stats_config_dict, stats_dir)

        # Create organvolumes section, if there are inverted labels
        if all(os.path.isdir(x) for x in [mut_inverted_labels, wt_inverted_labels]):
            org_config = {}   #'organvolumes'
            org_config['wt'] = wt_inverted_labels
            org_config['mut'] = mut_inverted_labels
            stats_config_dict['data']['organvolumes'] = org_config

        with open(stats_meta_path, 'w') as fh:
            fh.write(yaml.dump(stats_config_dict.to_dict(), default_flow_style=False))
        return stats_meta_path

    def add_intensity_stats_config(self, stats_config_dict, stats_dir):
        int_config = addict.Dict()  # intensity
        # If littermate controls are included with the mutants, we need to create a subset list of mutants to use

        wt_int = self.get_last_reg_stage(self.wt_path_maker)
        wt_intensity_dir = relpath(wt_int, stats_dir)
        mut_int = self.get_last_reg_stage(self.mut_path_maker)
        mut_intensity_dir = relpath(mut_int, stats_dir)

        intensity_normalisation_roi = self.wt_config.get('normalisation_roi')

        if not os.path.exists(wt_int):
            raise ValueError('cannot find wt intensity data\n{}'.format(wt_int))

        int_config['wt'] = wt_intensity_dir
        int_config['mut'] = mut_intensity_dir
        int_config['normalisation_roi'] = intensity_normalisation_roi
        stats_config_dict['data']['intensity'] = int_config

    def add_deformations_stats_config(self, stats_config_dict, stats_dir):
        # Jacobians and deformations can be generated at differnt scales eg 192-8, or 64-8
        # For now just do stats for the first scale in the list. Can craft a stats.yaml by hand if multiple scales
        # should be analysed

        def_config_entry = self.mut_config.get('generate_deformation_fields')
        if def_config_entry:
            first_def_scale = def_config_entry.keys()[0]
        wt_deformation_dir = relpath(join(self.wt_path_maker.get(DEFORMATION_DIR)), stats_dir)
        mut_deformation_dir = relpath(join(self.mut_path_maker.get(DEFORMATION_DIR)), stats_dir)
        wt_jacobian_dir = relpath(join(self.wt_path_maker.get(JACOBIAN_DIR), first_def_scale), stats_dir)
        mut_jacobian_dir = relpath(join(self.mut_path_maker.get(JACOBIAN_DIR), first_def_scale), stats_dir)

        if self.wt_config.get('generate_deformation_fields'):
            for deformation_id in self.wt_config['generate_deformation_fields'].keys():
                deformation_id = str(deformation_id)  # yaml interperets numbers with underscores as ints
                wt_deformation_scale_dir = join(wt_deformation_dir, deformation_id)
                wt_jacobian_scale_dir = wt_jacobian_dir #join(wt_jacobian_dir, deformation_id)
                mut_deformation_scale_dir = join(mut_deformation_dir, deformation_id)
                mut_jacobian_scale_dir = mut_jacobian_dir #join(mut_jacobian_dir, deformation_id)

                jacobians_scale_config = {
                    'wt': wt_jacobian_scale_dir,
                    'mut': mut_jacobian_scale_dir,
                }
                #returnFor now don't do dfeormations as it breaks in lama_stats
                # deformations_scale_config = {
                #     'wt_dir': wt_deformation_scale_dir,
                #     'mut_dir': mut_deformation_scale_dir,
                # }

                #
                # stats_config_dict['data']['deformations_' + deformation_id] = deformations_scale_config
                stats_config_dict['data']['jacobians_' + deformation_id] = jacobians_scale_config
                # Create a config file for the stats module to use

    def run_registration(self, config):
        lama.RegistraionPipeline(config, create_modified_config=False)

    def get_config(self, wt_config_path, mut_in_dir):
        """
        Gets the config file that was used for the wildtype registration.
        Copies it and fills out the mutant-specific entries eg relative paths to wt target, mask etc from the mut config
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

        def add_new_relative_path_to_mutant_config(config_parameter):
            if not wt_config.get(config_parameter):
                return
            wt_parameter_path = join(wt_config_dir, wt_config[config_parameter])
            parameter_path_rel_to_mut_config = relpath(wt_parameter_path, mut_config_dir)
            mutant_config[config_parameter] = parameter_path_rel_to_mut_config

        map(add_new_relative_path_to_mutant_config,
            ['label_map', 'label_names', 'isosurface_dir', 'fixed_volume', 'fixed_mask', 'staging_volume'])

        mutant_config['pad_dims'] = wt_config['pad_dims']

        return wt_config, wt_config_dir, mutant_config, mut_config_dir

    def get_last_reg_stage(self, path_object):
        last_reg_id = self.wt_config['registration_stage_params'][-1]['stage_id']
        reg_path = join(path_object.get('root_reg_dir'), last_reg_id)

        return reg_path


if __name__ == '__main__':

    # Log all uncaught exceptions
    sys.excepthook = common.excepthook_overide

    import argparse

    parser = argparse.ArgumentParser("MRC Harwell registration pipeline")
    parser.add_argument('-c', '--config', dest='wt_config', help='Config file of the wildtype run (YAML format)', required=True)
    parser.add_argument('-i', '--input', dest='in_dir', help='directory containing input volumes', required=False)
    parser.add_argument('-p', '--proj-dir', dest='mut_proj_dir', help='directory to put results', required=True)
    parser.add_argument('-n1', '--specimen_n=1', dest='n1', help='Do one mutant against many wts analysis?', default=False)
    parser.add_argument('-w', '--wildtpe_list', dest='wt_list', help='List of volume names that defines a subset of wt volumes to use', default=False)
    parser.add_argument('-l', '--littermate_bsaelines', dest='littermates', help='csv file with a list of liitermate filnames', default=False)
    parser.add_argument('-a', '--autostage', dest='autostage', help='use -a if baselines are to be chosen automatically. If -a not used, all baselines will be used', default=True)
    args, _ = parser.parse_known_args()
    PhenoDetect(args.wt_config, args.mut_proj_dir, args.wt_list, args.autostage, args.littermates, args.in_dir)

    args = parser.parse_args()

