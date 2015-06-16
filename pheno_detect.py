#!/usr/bin/python

"""
**************
phenodetect.py
**************


point towards a previous registration output of wildtypes. extract all the relevant settings and perform registration
of mutants.
"""
import os
import copy
import logging

import yaml

import mrch_regpipeline
from reg_stats import reg_stats

LOG_FILE = 'phenotype_detection.log'
LOG_MODE = logging.DEBUG

class PhenoDetect(object):
    def __init__(self, wt_config_path, in_dir, out_dir):

        logfile = os.path.join(out_dir, LOG_FILE)
        logging.basicConfig(filename=logfile, level=LOG_MODE, filemode="w")
        self.mutant_config, self.wt_output_metadata = self.get_config(wt_config_path, in_dir, out_dir)
        self.out_dir = self.mutant_config['output_dir']
        self.run_registration(self.mutant_config)

        self.stats_outdir = os.path.join(self.out_dir, 'stats')
        mrch_regpipeline.mkdir_force(self.stats_outdir)

        mutant_output_filename = os.path.join(self.out_dir, self.mutant_config['output_metadata_file'])
        self.mutant_output_metadata = yaml.load(open(mutant_output_filename, 'r'))

        self.intensity_analysis()
        self.jacobian_analysis()
        self.deformation_analysis()

    def intensity_analysis(self):
        out_file = os.path.join(self.stats_outdir, 'intensity.nrrd')
        wts = self.wt_output_metadata.get('normalized_registered')
        mutants = self.mutant_output_metadata['normalized_registered']
        if not wts:
            logging.warn("intensity analysis not performed as no normalised registered images were found")
            return
        reg_stats(wts, mutants, 'int', out_file, self.wt_output_metadata['fixed_mask'])

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
        mrch_regpipeline.RegistraionPipeline(config)

    def get_config(self, wt_config_path, in_dir, out_dir):
        wt_config = yaml.load(open(wt_config_path, 'r'))
        mutant_config = copy.deepcopy(wt_config)

        mutant_config['output_dir'] = os.path.abspath(out_dir)
        wt_out = wt_config['output_dir']
        wt_metadata_filename = wt_config['output_metadata_file']
        wt_output_metadata = yaml.load(open(os.path.join(wt_out, wt_metadata_filename), 'r'))

        mutant_config['inputvolumes_dir'] = os.path.abspath(in_dir)
        mutant_config['fixed_volume'] = wt_output_metadata['fixed_volume']
        mutant_config['fixed_mask'] = wt_output_metadata.get('fixed_mask')   # not required. will be set to None if not present

        return mutant_config, wt_output_metadata




if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("MRC Harwell phenotype detection tool")
    parser.add_argument('-w', dest='wt_config', help='Config file of the wildtype run (YAML format)', required=True)
    parser.add_argument('-i', dest='in_dir', help='directory containing input volumes', required=True)
    parser.add_argument('-o', dest='out_dir', help='directory to put results', required=True)
    args = parser.parse_args()
    PhenoDetect(args.wt_config, args.in_dir, args.out_dir)
