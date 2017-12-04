import nose
import os
from os.path import join, realpath, dirname, isdir
import sys
sys.path.insert(0, join(dirname(__file__), '..'))
import logging
import run_lama_stats

"""
Run all the config files in the test config directory 
Currently just running and making sure there's no uncaught exceptions.
TODO: check that the correct output is generated too
"""


current_dir = dirname(realpath(__file__))

CONFIG_DIR = join(current_dir, 'test_data', 'config_files')

def test_all():
	config = join(CONFIG_DIR, 'all_specimens.yaml')
	run_lama_stats.run(config)