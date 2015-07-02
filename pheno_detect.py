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
import copy
import logging

import yaml

import mrch_regpipeline
from reg_stats import reg_stats
import common

MUTANT_CONFIG = 'mutant_config.yaml'
"""str: Location to save the genrated config file for registering the mutants"""

LOG_FILE = 'phenotype_detection.log'

INTENSITY_DIR = 'normalized_registered'
"""str: directory to save the normalised registered images to"""

JACOBIAN_DIR = 'jacobians'
"""str: directory to save the determinant of the jacobians to"""

DEFORMATION_DIR = 'deformation_fields'
"""str: directory to save the deformation fileds to"""

class PhenoDetect(object):
    """Phenotype detection .
    """
    def __init__(self, wt_config_path, proj_dir, in_dir):

        self.proj_dir = os.path.abspath(proj_dir)

        logfile = os.path.join(self.proj_dir, LOG_FILE)
        common.init_log(logfile, "Phenotype detectoin pipeline")

        self.wt_config_dir = os.path.split(os.path.abspath(wt_config_path))[0]

        self.mutant_config, self.wt_output_metadata = self.get_config(wt_config_path, in_dir)
        self.out_dir = os.path.join(self.proj_dir, self.mutant_config['output_dir'])
        self.mutant_config['output_dir'] = self.out_dir
        mutant_config_path = self.write_config()

        if not self.wt_output_metadata.get('fixed_mask'):
            self.fixed_mask = None
            logging.warn('WT fixed mask not present. Optimal results will not be obtained')
        else:
            self.fixed_mask = os.path.join(self.proj_dir, self.wt_output_metadata['fixed_mask'])

        self.run_registration(mutant_config_path)

        self.stats_outdir = os.path.join(self.out_dir, 'stats')
        mrch_regpipeline.mkdir_force(self.stats_outdir)

        mutant_output_filename = os.path.join(self.out_dir, self.mutant_config['output_metadata_file'])
        self.mutant_output_metadata = yaml.load(open(mutant_output_filename, 'r'))

        common.log_time('Stats analysis started')
        self.run_analysis('intensity', INTENSITY_DIR, 'scalar')
        self.run_analysis('jacobians', JACOBIAN_DIR, 'scalar')
        self.run_analysis('deformations', DEFORMATION_DIR, 'vector')

    def write_config(self):
        config_path = os.path.join(self.proj_dir, MUTANT_CONFIG)
        with open(config_path, 'w') as fh:
            fh.write(yaml.dump(self.mutant_config, default_flow_style=False))
        return config_path

    def run_analysis(self, out_name, directory, data_type):
        out_dir = os.path.join(self.stats_outdir, out_name)
        try:
            common.mkdir_force(out_dir)
        except OSError as e:
            logging.warn("Can create directory: {}".format(e))
            return

        wts = self.wt_output_metadata.get(directory)
        if not wts:
            logging.warn("Can't find path '{}' so can't create {}".format(directory, out_name))
            return
        wt_abs = os.path.join(self.wt_config_dir, wts)
        mutants = self.mutant_output_metadata[directory]
        mut_abs = os.path.join(self.proj_dir, mutants)
        reg_stats(wt_abs, mut_abs, data_type, out_dir, self.fixed_mask, n_eq_1=True)

    def jacobian_analysis(self):
        out_file = os.path.join(self.stats_outdir, 'jacobians.nrrd')
        wts = self.wt_output_metadata.get('jacobians')
        mutants = self.mutant_output_metadata['jacobians']
        if not wts:
            logging.warn("jacobian analysis not performed")
            return
        reg_stats(wts, mutants, 'jac', out_file, self.wt_output_metadata['fixed_mask'])

    def deformation_analysis(self):
        out_file = os.path.join(self.stats_outdir, 'deformations.nrrd')
        wts = self.wt_output_metadata.get('deformation_fields')
        mutants = self.mutant_output_metadata['deformation_fields']
        if not wts:
            logging.warn("deformation analysis not performed")
            return
        reg_stats(wts, mutants, 'def', out_file, self.wt_output_metadata['fixed_mask'])

    def run_registration(self, config):
        mrch_regpipeline.RegistraionPipeline(config, phenotyping=True)

    def get_config(self, wt_config_path, in_dir):
        wt_config = yaml.load(open(wt_config_path, 'r'))
        mutant_config = copy.deepcopy(wt_config)

        wt_metadata_filename = os.path.join(self.wt_config_dir, 'out',  wt_config['output_metadata_file'])
        wt_output_metadata = yaml.load(open(wt_metadata_filename, 'r'))

        mutant_config['inputvolumes_dir'] = in_dir
        mutant_config['fixed_volume'] = wt_output_metadata['fixed_volume']
        mutant_config['fixed_mask'] = wt_output_metadata.get('fixed_mask')   # not required. will be set to None if not present
        mutant_config['pad_dims'] = wt_output_metadata.get('pad_dims')
        mutant_config['wt_proj_dir'] = os.path.dirname(wt_config_path)

        return mutant_config, wt_output_metadata


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("MRC Harwell phenotype detection tool")
    parser.add_argument('-c', '--config', dest='wt_config', help='Config file of the wildtype run (YAML format)', required=True)
    parser.add_argument('-i', '--input', dest='in_dir', help='directory containing input volumes', required=True)
    parser.add_argument('-p', '--proj-dir', dest='proj_dir', help='directory to put results', required=True)
    #parser.add_argument('-l', '--labelmap', dest='proj_dir', help='directory to put results', default=None)
    # parser.add_argument('-l', '--labelmap', dest='proj_dir', help='directory to put results', default=None)

    args = parser.parse_args()
    PhenoDetect(args.wt_config, args.proj_dir, args.in_dir)
