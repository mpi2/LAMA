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
from nose.tools import assert_raises, eq_, nottest
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


# @nottest
def test_p_thresholds():
    """
    Testing the p_thresholds calculation
    -------
    TODO: Add more tests for different cases
    """
    import copy
    import pandas as pd

    # These simulate null distributions from 100 baselines for two organs
    null = pd.DataFrame.from_dict({
        '1': np.linspace(0.01, 1, 100),
        '2': np.linspace(0.01, 1, 100)})

    # These simulate alternative distributions for two organs from 40 lines
    alt = pd.DataFrame.from_dict({
        '1': np.linspace(0.00000001, 0.1, 40),
        '2': np.linspace(0.9, 1, 40)})

    thresh = p_thresholds.get_thresholds(null, alt)

    eq_(thresh.loc[1, 'p_thresh'],  0.02)  # Gives a p-value threshold of 0.02
    eq_(thresh.loc[2, 'p_thresh'], 1.0)    # Gives a p-value threshold of 1.0 as there are no low p-values in the alt distribution




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

