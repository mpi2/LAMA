"""
Test the permutation-based stats pipeline

Usage
-----

To run just this test:
    $ cd lama/
    $ nosetests tests/_test_run_permutation_stats.py
"""


import sys
from pathlib import Path
import numpy as np
from pathlib import Path
from nose.tools import assert_raises, nottest
from lama.stats.permutation_stats import run_permutation_stats
from lama.stats.permutation_stats import p_thresholds
import pandas as pd
# Import from __init__
from tests import test_data_root, registration_data_root


outdir = test_data_root / 'test_output'

# @nottest
def test_permutation_stats():
    """
    Run the whole permutation based stats pipeline.
    Currently this just checks to see if the pipeline completes without any errors.
    """
    num_perms = 20  # Would do 1000 or more normally
    wt_dir = registration_data_root / 'baseline'
    mut_dir = registration_data_root / 'mutant'
    perm_out_dir = outdir / 'organ_vols_permutation' # Intermediate results go here. Permutation distributions etc.

    run_permutation_stats.run(wt_dir, mut_dir, perm_out_dir, num_perms, log_dependent=False)


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
    alt_file = Path('/home/neil/git/lama/tests/test_data/stats_test_data/test_output/organ_vols_permutation/alt_line_dist_pvalues.csv')
    thresholds_file = Path('/home/neil/git/lama/tests/test_data/stats_test_data/test_output/organ_vols_permutation/line_organ_p_thresholds.csv')
    mutant_dir = Path('/home/neil/git/lama/tests/test_data/registration_test_data/mutant')

    thresholds = pd.read_csv(thresholds_file, index_col=0)
    alt = pd.read_csv(alt_file, index_col=0)

    run_permutation_stats.annotate(thresholds, alt, mutant_dir)

    # # Specimens
    alt_file = Path('/home/neil/git/lama/tests/test_data/stats_test_data/test_output/organ_vols_permutation/alt_specimen_dist_pvalues.csv')
    thresholds_file = Path('/home/neil/git/lama/tests/test_data/stats_test_data/test_output/organ_vols_permutation/specimen_organ_p_thresholds.csv')
    mutant_dir = Path('/home/neil/git/lama/tests/test_data/registration_test_data/mutant')

    thresholds = pd.read_csv(thresholds_file, index_col=0)
    alt = pd.read_csv(alt_file, index_col=0)

    run_permutation_stats.annotate(thresholds, alt, mutant_dir)

