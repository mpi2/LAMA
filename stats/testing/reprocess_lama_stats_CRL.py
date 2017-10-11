#!/usr/bin/python

"""
11/10/2017

Reruning the current lines with lama_stats.py bu using CRl as a fixed effect
lama_stats expects a config file
"""

from os.path import expanduser, join, dirname
home_dir = expanduser('~')
import sys
import yaml
sys.path.insert(0, join(dirname(__file__), '..'))
import lama_stats as stats
import common

config_dict = {'data':             {'jacobians_192_to_12': {'mut_dir': '../jacobians/192_to_12',
                                                            'wt_dir':  '../../../../output/jacobians/192_to_12'}},
               'fixed_mask':       '../../../../output/padded_target/mask_cleaned.nrrd',
               'formulas':         ['voxel ~ genotype + crl'],
               'label_map_path':   '../../../../output/padded_target/seg_compressed.nrrd',
               'mut_staging_file': '../staging_info.csv',
               'n1':               True,
               'voxel_size':       14.0,
               'wt_staging_file':  '../../../../output/staging_info.csv'}

lines_list_path = join(home_dir, 'bit/LAMA_results/E14.5/paper_runs/mutant_runs/reprocess_lama_stats_crl_fe.csv')
root_dir = join(home_dir, 'bit/LAMA_results/E14.5/paper_runs/mutant_runs')


def run_stats(line_name):
	config_dict['project_name'] = line_name
	outdir = join(root_dir, line_name, 'output', 'stats_crl')
	try:
		common.mkdir_force(outdir)
		config_path = join(outdir, 'stats_jac_crl.yaml')
	except (IOError, OSError) as e:
		print('{} did not work has it got an output folder?\n{}'.format(line_name, str(e)))

	with open(config_path, 'w') as fh:
		fh.write(yaml.dump(config_dict))
	try:
		stats.run(config_path)
	except Exception as e:
		print('{} failed\n{}.'.format(line_name, e))


lines = []
with open(lines_list_path, 'r') as fh:
	for line in fh:
		lines.append(line.strip())

for line in lines:
	run_stats(line)
