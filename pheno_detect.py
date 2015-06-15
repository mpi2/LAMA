#!/usr/bin/python

"""
phenodetect.py

point towards a previous registration output of wildtypes. extract all the relevant settings and perform registration
of mutants.
"""

import yaml
import copy
import os
from mrch_regpipeline import RegistraionPipeline

WT_OUT_PATHS = 'out_paths.yaml'

class PhenoDetect(object):
    def __init__(self, wt_config_path, in_dir, out_dir):

        self.mutant_config = self.get_config(wt_config_path, in_dir, out_dir)
        self.run_registration(self.mutant_config)

    def run_registration(self, config):
        RegistraionPipeline(config)

    def get_config(self, wt_config_path, in_dir, out_dir):
        wt_config = yaml.load(open(wt_config_path, 'r'))
        mutant_config = copy.deepcopy(wt_config)

        mutant_config['output_dir'] = os.path.abspath(out_dir)
        wt_out = wt_config['output_dir']
        wt_paths = yaml.load(open(os.path.join(wt_out, WT_OUT_PATHS), 'r'))

        mutant_config['inputvolumes_dir'] = os.path.abspath(in_dir)
        mutant_config['fixed_volume'] = wt_paths['fixed_volume']
        mutant_config['fixed_mask'] = wt_paths.get('fixed_mask')   # not required. will be set to None if not present

        return mutant_config




if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("MRC Harwell phenotype detection tool")
    parser.add_argument('-w', dest='wt_config', help='Config file of the wildtype run (YAML format)', required=True)
    parser.add_argument('-i', dest='in_dir', help='directory containing input volumes', required=True)
    parser.add_argument('-o', dest='out_dir', help='directory to put results', required=True)
    args = parser.parse_args()
    PhenoDetect(args.wt_config, args.in_dir, args.out_dir)
