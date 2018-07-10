#!/usr/bin/python

"""
11/10/2017

Reruning the current lines with run_lama_stats.py bu using CRl as a fixed effect
lama_stats expects a config file
"""

from os.path import expanduser, join, dirname, isdir

home_dir = expanduser('~')
import sys
import yaml

sys.path.insert(0, join(dirname(__file__), '..'))
import run_lama_stats as stats
import common
import shutil
config_dict = {'data': {'organvolumes': {'mut': '../inverted_labels/similarity',
                                                'wt': '../../../../../output/inverted_labels/similarity'}},
               'fixed_mask': '../../../../../output/padded_target/mask_june_18.nrrd',
               'formulas': ['voxel ~ genotype + crl'],
               'label_map': '../../../../../output/padded_target/v24_dev30.nrrd',
               'label_names':  '../../../../../output/padded_target/dev24_atlas_terms_with_organ_system.csv',
               'mut_staging_file': '../staging_info.csv',
               'n1': False,
               'voxel_size': 14.0,
               'wt_staging_file': '../../../../../output/staging_info.csv',
               'littermate_pattern': '_wt_',
               'use_auto_staging': False}

lines_list_path = join(home_dir, 'bit/LAMA_results/E14.5/paper_runs/mutant_runs/280618_analysed_lines/lines.csv')
root_dir = join(home_dir, 'bit/LAMA_results/E14.5/paper_runs/mutant_runs/280618_analysed_lines')
log_path = join(home_dir, 'bit/LAMA_results/E14.5/paper_runs/mutant_runs/090718_analysed_lines.log')


def run_stats(line_name):
    with open(log_path, 'a') as logger:
        config_dict['project_name'] = line_name

        outdir = join(root_dir, line_name, 'output', 'stats')


        if isdir(outdir):
            shutil.rmtree(outdir)
        common.mkdir_force(outdir)

        config_path = join(outdir, 'stats_organ_crl.yaml')
        with open(config_path, 'w') as fh:
            fh.write(yaml.dump(config_dict))
        try:
            stats.run(config_path)
        except Exception as e:
            msg = '{} failed\n{}.'.format(line_name, e)
            print(msg)
            logger.write(msg)




lines = []
with open(lines_list_path, 'r') as fh:
    for line in fh:
        lines.append(line.strip())

for line in lines:
    print('doing {}'.format(line))
    run_stats(line)
