import nose
from nose.tools import nottest, raises
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

To run these tests, the test data needs to be fected from bit/dev/lama_stats_test_data
In future we should put this in a web-accessible place
"""


current_dir = dirname(realpath(__file__))

CONFIG_DIR = join(current_dir, 'test_data', 'config_files')


def test_all():
	config = join(CONFIG_DIR, 'all_specimens.yaml')
	run_lama_stats.run(config)


def test_all():
	config = join(CONFIG_DIR, 'all_specimens.yaml')
	run_lama_stats.run(config)


@raises(IOError)
def test_missing_config():
	config = join(CONFIG_DIR, 'all_specimens.flaml')
	run_lama_stats.run(config)