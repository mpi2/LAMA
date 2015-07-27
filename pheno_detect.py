#!/usr/bin/env python


"""pheno_detect.py

Phenotype detection component of the MRC Harwell registration pipeline

Example
-------

    $ pheno_detect.py -c wt_registration_config_file_path -i dir_with_input_volumes -p project_directory_path

Notes
-----

TODO: add label maps for inversion

"""

import os
from os.path import join, relpath
import copy
import logging

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
    """Phenotype detection

    TODO
    ----
    This does not really need to be a class
    """
    def __init__(self, wt_config_path, proj_dir, in_dir, n1=True):

        self.proj_dir = os.path.abspath(proj_dir)
        self.n1 = n1

        logfile = join(self.proj_dir, LOG_FILE)
        common.init_log(logfile, "Phenotype detectoin pipeline")

        self.wt_config_dir = os.path.split(os.path.abspath(wt_config_path))[0]
        self.wt_output_metadata_dir = ''  # Dir containing the metadat file of the wt run. Gets set in get_config (change)
        #self.mut_output_metadata_dir = ''

        self.mutant_config, self.wt_output_metadata = self.get_config(wt_config_path, in_dir)
        self.out_dir = join(self.proj_dir, self.mutant_config['output_dir'])
        self.mutant_config['output_dir'] = self.out_dir
        mutant_config_path = self.write_config()

        if not self.wt_output_metadata.get('fixed_mask'):
            self.fixed_mask = None
            logging.warn('WT fixed mask not present. Optimal results will not be obtained')
        else:
            self.fixed_mask = join(self.wt_output_metadata_dir, self.wt_output_metadata['fixed_mask'])

        self.run_registration(mutant_config_path)

        mutant_output_filename = join(self.out_dir, self.mutant_config['output_metadata_file'])
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
        After getting the wildtype registration config and substituting mutnat-specific info, write out a mutant
        registration config file
        """
        config_path = join(self.proj_dir, MUTANT_CONFIG)
        with open(config_path, 'w') as fh:
            fh.write(yaml.dump(self.mutant_config, default_flow_style=False))
        return config_path

    def write_stats_config(self):
        """
        Writes a yaml config file for use by the reg_stats.py module to use. Provides paths to data and some options
        """

        wt_intensity_dir = relpath(join(self.wt_output_metadata_dir, self.wt_output_metadata.get(INTENSITY_DIR)), self.out_dir)
        wt_deformation_dir = relpath(join(self.wt_output_metadata_dir, self.wt_output_metadata.get(DEFORMATION_DIR)), self.out_dir)
        wt_jacobian_dir = relpath(join(self.wt_output_metadata_dir, self.wt_output_metadata.get(JACOBIAN_DIR)), self.out_dir)

        mut_intensity_dir = relpath(join(self.out_dir, self.mutant_output_metadata[INTENSITY_DIR]), self.out_dir)
        mut_deformation_dir = relpath(join(self.out_dir, self.mutant_output_metadata[DEFORMATION_DIR]), self.out_dir)
        mut_jacobian_dir = relpath(join(self.out_dir, self.mutant_output_metadata[JACOBIAN_DIR]), self.out_dir)

        fixed_mask = relpath(self.fixed_mask, self.out_dir)

        stats_meta_path = join(self.out_dir, STATS_METADATA_PATH)

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
        mrch_regpipeline.RegistraionPipeline(config, phenotyping=True)

    def get_config(self, wt_config_path, in_dir):
        wt_config = yaml.load(open(wt_config_path, 'r'))
        mutant_config = copy.deepcopy(wt_config)

        wt_metadata_filename = join(self.wt_config_dir, 'out',  wt_config['output_metadata_file'])
        self.wt_output_metadata_dir = os.path.dirname(wt_metadata_filename)

        wt_output_metadata = yaml.load(open(wt_metadata_filename, 'r'))

        mutant_config['inputvolumes_dir'] = in_dir
        mutant_config['fixed_volume'] = join(self.wt_output_metadata_dir, wt_output_metadata['fixed_volume'])
        mutant_config['fixed_mask'] = join(self.wt_output_metadata_dir, wt_output_metadata.get('fixed_mask'))   # not required. will be set to None if not present
        mutant_config['pad_dims'] = wt_output_metadata.get('pad_dims')
        mutant_config['wt_proj_dir'] = os.path.dirname(wt_config_path)

        return mutant_config, wt_output_metadata


if __name__ == '__main__':

    import argparse
    import sys

    parser = argparse.ArgumentParser("MRC Harwell registration pipeline")
    parser.add_argument('-c', '--config', dest='wt_config', help='Config file of the wildtype run (YAML format)', required=True)
    parser.add_argument('-i', '--input', dest='in_dir', help='directory containing input volumes', required=True)
    parser.add_argument('-p', '--proj-dir', dest='proj_dir', help='directory to put results', required=True)
    parser.add_argument('-n1', '--specimen_n=1', dest='n1', help='Do one mutant against many wts analysis?', default=False)
    args, _ = parser.parse_known_args()
    PhenoDetect(args.wt_config, args.proj_dir, args.in_dir)




    # parser.add_argument('-l', '--labelmap', dest='proj_dir', help='directory to put results', default=None)

    args = parser.parse_args()




