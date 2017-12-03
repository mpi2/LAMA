import nose
import os
from os.path import join, realpath, dirname, isdir
import sys
sys.path.insert(0, join(dirname(__file__), '..'))
import logging
import run_lama_stats

"""
Given a set of lama_stats config files run and check for errors and the desired output being present
"""


current_dir = dirname(realpath(__file__))

CONFIG_DIR = join(current_dir, 'test_data', 'config_files')

config_files =['config_1.yaml'
]

def test_stats():
	for c in config_files:
		# Skip output directories
		config = join(CONFIG_DIR, c)
		run_lama_stats.run(config)