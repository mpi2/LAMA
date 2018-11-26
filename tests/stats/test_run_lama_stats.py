from os.path import join, dirname
import sys
sys.path.insert(0, join(dirname(__file__), '../..', 'stats'))
sys.path.insert(0, join(dirname(__file__), '../..', 'stats', 'newstats'))
from logzero import logger as logging
import os
import numpy as np
from pathlib import Path
from nose.tools import assert_raises, nottest
import run_lama_stats
import run_permutation_stats
import p_thresholds
from . import CONFIG_DIR

import pandas as pd

"""
Run all the config files in the test config directory 
Currently just running and making sure there's no uncaught exceptions.
TODO: check that the correct output is generated too

To run these tests, the test data needs to be fechted from bit/dev/lama_stats_test_data
In future we should put this in a web-accessible place
"""
@nottest
def test_organ_vols():
    config = join(CONFIG_DIR, 'organ_vols.yaml')
    run_lama_stats.run(config)

@nottest
def test_all():
    config = join(CONFIG_DIR, 'all_specimens.yaml')
    run_lama_stats.run(config)


@nottest
def test_glcm():
    config = join(CONFIG_DIR, 'test_glcm.yaml')
    run_lama_stats.run(config)

@nottest
def test_from_arbitrary_directory():
    config = join(CONFIG_DIR, 'all_specimens.yaml')
    os.chdir(os.path.expanduser('~'))
    run_lama_stats.run(config)

@nottest
def test_missing_config():
    config = join(CONFIG_DIR, 'all_specimens.flaml')
    assert_raises(Exception, run_lama_stats.run, config)

@nottest
def test_misc():
    config = join(CONFIG_DIR, 'misc_test.yaml')
    run_lama_stats.run(config)


@nottest
def test_incorrect_staging_file():
    pass


# @nottest
def test_permutation_stats():
    """
    Run the whole permutation based stats pipeline.
    Just looking to see whether it completed without any errore
    """
    # TODO: remove hardoding of these paths

    num_perms = 10  # Would do 1000 or more normally
    wt_dir = Path('/home/neil/git/lama/tests/test_data/registration_test_data/baseline')
    mut_dir = Path('/home/neil/git/lama/tests/test_data/registration_test_data/mutant')
    out_dir = Path('/home/neil/git/lama/tests/test_data/stats_test_data/test_output/organ_vols_permutation')


    run_permutation_stats.run(wt_dir, mut_dir, out_dir, num_perms, log_dependent=True)


@nottest
def test_p_thresholds():
    """
    Testing the p_thresholds calculation
    -------
    TODO: Add more tests for different cases
    """
    import copy
    import pandas as pd
    # create a null distribution for two example organs with 10 permutation
    null_1 = np.linspace(0.01, 1, 100)  # array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1 , ..... 1.0
    null_2 = np.linspace(0.01, 1, 100)

    alt_1 = copy.copy(null_1) / 2
    alt_2 = copy.copy(null_2) / 10

    null = pd.DataFrame.from_dict({'1': null_1,
                                   '2': null_2})

    alt = pd.DataFrame.from_dict({'1': alt_1,
                                  '2': alt_2})

    thresh = p_thresholds.get_thresholds(null, alt)
    print(thresh)  # TODO: Check it's giving correct values back

@nottest
def test_annotate():
    # Lines
    # alt_file = Path('/home/neil/git/lama/tests/test_data/stats_test_data/test_output/organ_vols_permutation/alt_line_dist_pvalues.csv')
    # thresholds_file = Path('/home/neil/git/lama/tests/test_data/stats_test_data/test_output/organ_vols_permutation/line_organ_p_thresholds.csv')
    # mutant_dir = Path('/home/neil/git/lama/tests/test_data/registration_test_data/mutant')
    #
    # thresholds = pd.read_csv(thresholds_file, index_col=0)
    # alt = pd.read_csv(alt_file, index_col=0)
    #
    # run_permutation_stats.annotate(thresholds, alt, mutant_dir)

    # # Specimens
    alt_file = Path('/home/neil/git/lama/tests/test_data/stats_test_data/test_output/organ_vols_permutation/alt_specimen_dist_pvalues.csv')
    thresholds_file = Path('/home/neil/git/lama/tests/test_data/stats_test_data/test_output/organ_vols_permutation/specimen_organ_p_thresholds.csv')
    mutant_dir = Path('/home/neil/git/lama/tests/test_data/registration_test_data/mutant')

    thresholds = pd.read_csv(thresholds_file, index_col=0)
    alt = pd.read_csv(alt_file, index_col=0)

    run_permutation_stats.annotate(thresholds, alt, mutant_dir)

