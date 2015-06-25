#!/usr/bin/python

"""
**************
phenodetect.py
**************

#todo: send the pipe yaml file intead of dict
point towards a previous registration output of wildtypes. extract all the relevant settings and perform registration
of mutants.
"""
import os
import copy
import logging

import yaml

import mrch_regpipeline
from reg_stats import reg_stats

MUTANT_CONFIG = 'mutant_config.yaml'
LOG_FILE = 'phenotype_detection.log'
INTENSITY_DIR = 'normalized_registered'
JACOBIAN_DIR = 'jacobians'
DEFORMATION_DIR = 'deformation_fields'

LOG_MODE = logging.DEBUG

class PhenoDetect(object):
    def __init__(self, wt_config_path, proj_dir, in_dir):

        self.proj_dir = os.path.abspath(proj_dir)

        logfile = os.path.join(self.proj_dir, LOG_FILE)
        logging.basicConfig(filename=logfile, level=LOG_MODE, filemode="w")

        self.wt_config_dir = os.path.split(os.path.abspath(wt_config_path))[0]

        self.mutant_config, self.wt_output_metadata = self.get_config(wt_config_path, in_dir)
        self.out_dir = os.path.join(self.proj_dir, self.mutant_config['output_dir'])
        self.mutant_config['output_dir'] = self.out_dir
        mutant_config_path = self.write_config()

        self.fixed_mask = os.path.join(self.proj_dir, self.wt_output_metadata['fixed_mask'])

        self.run_registration(mutant_config_path)

        self.stats_outdir = os.path.join(self.out_dir, 'stats')
        mrch_regpipeline.mkdir_force(self.stats_outdir)

        mutant_output_filename = os.path.join(self.out_dir, self.mutant_config['output_metadata_file'])
        self.mutant_output_metadata = yaml.load(open(mutant_output_filename, 'r'))

        self.run_analysis('intensity.nrrd', INTENSITY_DIR, 'scalar')
        self.run_analysis('jacobians.nrrd', JACOBIAN_DIR, 'scalar')
        self.run_analysis('deformations.nrrd', DEFORMATION_DIR, 'vector')

    def write_config(self):
        config_path = os.path.join(self.proj_dir, MUTANT_CONFIG)
        with open(config_path, 'w') as fh:
            fh.write(yaml.dump(self.mutant_config, default_flow_style=False))
        return config_path

    def run_analysis(self, out_name, directory, data_type):
        out_file = os.path.join(self.stats_outdir, out_name)
        wts = self.wt_output_metadata.get(directory)
        wt_abs = os.path.join(self.proj_dir, wts)
        mutants = self.mutant_output_metadata[directory]
        mut_abs = os.path.join(self.proj_dir, mutants)

        if not wts:
            logging.warn("intensity analysis not performed as no normalised registered images were found")
            return
        reg_stats(wt_abs, mut_abs, data_type, out_file, self.fixed_mask)

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
    parser.add_argument('-c', dest='wt_config', help='Config file of the wildtype run (YAML format)', required=True)
    parser.add_argument('-i', dest='in_dir', help='directory containing input volumes', required=True)
    parser.add_argument('-p', dest='proj_dir', help='directory to put results', required=True)
    args = parser.parse_args()
    PhenoDetect(args.wt_config, args.proj_dir, args.in_dir)
