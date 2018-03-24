import nose
from nose.tools import nottest, raises
import os
from os.path import join, realpath, dirname, isdir
import sys

sys.path.insert(0, join(dirname(__file__), '..'))
from lama import RegistraionPipeline

"""
Run all the config files in the test config directory 
Currently just running and making sure there's no uncaught exceptions.
TODO: check that the correct output is generated too

To run these tests, the test data needs to be fechted from bit/dev/lama_stats_test_data
In future we should put this in a web-accessible place
"""

current_dir = dirname(realpath(__file__))

CONFIG_DIR = join(current_dir, 'test_data', 'config_files')


def test_all():
    config = join(CONFIG_DIR, 'all_specimens.yaml')
    run_lama_stats.run(config)
    # todo: read in the stats log and corm the min p q and min max t values are within a certain range
    # jacobians_config = join(current_dir, 'test_data', 'test_output')


def test_missing_config():
    config = join(CONFIG_DIR, 'all_specimens.flaml')
    try:
        run_lama_stats.run(config)
        assert False
    except IOError:
        assert True
