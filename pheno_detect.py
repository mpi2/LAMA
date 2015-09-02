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

import mrch_regpipeline
from reg_stats import reg_stats
from organvolume_stats import organvolume_stats
import common

VOLUME_CALCULATIONS_FILENAME = "organvolumes.csv"

MUTANT_CONFIG = 'mutant_config.yaml'
"""str: Location to save the genrated config file for registering the mutants"""

LOG_FILE = 'phenotype_detection.log'

INTENSITY_DIR = 'normalized_registered'
"""str: directory to save the normalised registered images to"""

JACOBIAN_DIR = 'jacobians'
"""str: directory to save the determinant of the jacobians to"""

DEFORMATION_DIR = 'deformation_fields'
"""str: directory to save the deformation fileds to"""

STATS_METADATA_HEADER = "This file can be run like: reg_stats.py -c stats.yaml"
STATS_METADATA_PATH = 'stats.yaml'


class PhenoDetect(object):
    """
    Phenotype detection

    """
    def __init__(self, wt_config_path, mut_proj_dir, in_dir, n1=True):
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

        self.in_dir = in_dir

        logfile = join(self.mut_proj_dir, LOG_FILE)
        common.init_log(logfile, "Phenotype detectoin pipeline")

        (self.wt_config, self.wt_config_dir,
         self.mut_config, self.mut_config_dir) = self.get_config(wt_config_path, in_dir)

        self.out_dir = join(self.mut_proj_dir, self.mut_config['output_dir'])
        self.mut_config['output_dir'] = self.out_dir


        self.mut_config_path = join(self.mut_proj_dir, MUTANT_CONFIG)
        self.write_config()

        self.run_registration(self.mut_config_path)

        mutant_output_filename = join(self.mut_proj_dir, self.out_dir, self.mutant_config['output_metadata_file'])
        self.mutant_output_metadata = yaml.load(open(mutant_output_filename, 'r'))

        common.log_time('Stats analysis started')

        stats_metadata_path = self.write_stats_config()
        reg_stats(stats_metadata_path)

        #Calculate organ volumes and do some stats on them

        wt_organ_vols = join(self.wt_output_metadata_dir,
                             self.wt_output_metadata['label_inversion_dir'],
                             VOLUME_CALCULATIONS_FILENAME)

        mut_organ_volumes = join(self.out_dir,
                                 self.wt_output_metadata['label_inversion_dir'],
                                 VOLUME_CALCULATIONS_FILENAME)
        organ_vol_stats_out = join(self.out_dir, 'organ_volume_stats.csv')

        organvolume_stats(wt_organ_vols, mut_organ_volumes, organ_vol_stats_out)

    def write_config(self):
        """
        After getting the wildtype registration config and substituting mutant-specific info, write out a mutant
        registration config file
        """

        replacements = {'inputvolumes_dir': self.mut_config['inputvolumes_dir'],
                        'fixed_volume': self.mut_config['fixed_volume'],
                        'fixed_mask': self.mut_config['fixed_mask'],
                        'label_map_path': self.mut_config['label_map_path'],
                        'pad_dims': self.mut_config['pad_dims'],
                        'organ_names': self.mut_config['organ_names']}

        mrch_regpipeline.replace_config_lines(self.mut_config_path, replacements)

    def write_stats_config(self):
        """
        Writes a yaml config file for use by the reg_stats.py module. Provides paths to data and some options
        """

        wt_intensity_dir = relpath(join(self.wt_output_metadata_dir, self.wt_output_metadata.get(INTENSITY_DIR)), self.out_dir)
        wt_deformation_dir = relpath(join(self.wt_output_metadata_dir, self.wt_output_metadata.get(DEFORMATION_DIR)), self.out_dir)
        wt_jacobian_dir = relpath(join(self.wt_output_metadata_dir, self.wt_output_metadata.get(JACOBIAN_DIR)), self.out_dir)

        mut_intensity_dir = relpath(join(self.out_dir, self.mutant_output_metadata[INTENSITY_DIR]), self.out_dir)
        mut_deformation_dir = relpath(join(self.out_dir, self.mutant_output_metadata[DEFORMATION_DIR]), self.out_dir)
        mut_jacobian_dir = relpath(join(self.out_dir, self.mutant_output_metadata[JACOBIAN_DIR]), self.out_dir)

        fixed_mask = relpath(self.fixed_mask, self.out_dir)

        stats_meta_path = join(self.mut_proj_dir, self.out_dir, STATS_METADATA_PATH)

        inverted_tform_dir = self.mutant_output_metadata['inverted_elx_dir']
        inverted_tform_config = join(inverted_tform_dir, "invert.yaml")

        # Create a metadat file for the stats module to use
        stats_metadata = {
            'fixed_mask': fixed_mask,
            'n1': self.n1,
            'data': {
                'registered_normalised':
                    {'datatype': 'scalar',
                     'wt': wt_intensity_dir,
                     'mut': mut_intensity_dir
                     },
                'glcm':
                    {'datatype': 'scalar',
                     'wt': wt_intensity_dir,
                     'mut': mut_intensity_dir
                     },
                'deformations':
                    {'datatype': 'vector',
                     'wt': wt_deformation_dir,
                     'mut': mut_deformation_dir
                     },
                'jacobians':
                    {'datatype': 'scalar',
                     'wt': wt_jacobian_dir,
                     'mut': mut_jacobian_dir
                     }
            },
            'inverted_tform_config': inverted_tform_config
        }

        with open(stats_meta_path, 'w') as fh:
            fh.write(yaml.dump(stats_metadata, default_flow_style=False))
        return stats_meta_path

    def run_registration(self, config):
        mrch_regpipeline.RegistraionPipeline(config)

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

        fixed_labelmap = join(wt_config_dir, wt_config['label_map_path'])
        fixed_labelmap_rel = relpath(fixed_labelmap, mut_config_dir)

        fixed_organ_names = join(wt_config_dir, wt_config['organ_names'])
        organ_names_rel = relpath(fixed_organ_names, mut_config_dir)

        mutant_config['fixed_volume'] = fixed_vol_rel
        mutant_config['fixed_mask'] = fixed_mask_rel
        mutant_config['pad_dims'] = wt_config['pad_dims']
        mutant_config['label_map_path'] = fixed_labelmap_rel
        mutant_config['organ_names'] = organ_names_rel

        return wt_config, wt_config_dir, mutant_config, mut_config_dir


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser("MRC Harwell registration pipeline")
    parser.add_argument('-c', '--config', dest='wt_config', help='Config file of the wildtype run (YAML format)', required=True)
    parser.add_argument('-i', '--input', dest='in_dir', help='directory containing input volumes', required=True)
    parser.add_argument('-p', '--proj-dir', dest='mut_proj_dir', help='directory to put results', required=True)
    parser.add_argument('-n1', '--specimen_n=1', dest='n1', help='Do one mutant against many wts analysis?', default=False)
    args, _ = parser.parse_known_args()
    PhenoDetect(args.wt_config, args.mut_proj_dir, args.in_dir)


    args = parser.parse_args()




