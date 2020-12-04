"""
Test the permutation-based stats pipeline

Usage
-----

These functions test the lama registration pipeline permutation stats module

Usage:  pytest -v -m "not notest" test_run_permutation_stats.py
"""


import sys
from pathlib import Path
import numpy as np
from pathlib import Path
import pytest

from nose.tools import assert_raises, eq_, ok_, nottest
from lama.stats.permutation_stats import run_permutation_stats
from lama.stats.permutation_stats import p_thresholds
import pandas as pd
# Import from __init__
from lama.tests.archive import test_data_root, registration_root, wt_registration_dir, mut_registration_dir, target_dir
from lama.common import ORGAN_VOLUME_CSV_FILE, STAGING_INFO_FILENAME

outdir = test_data_root / 'test_output'


def test_prepare_data(n_baselines=100, n_mutant_lines=20, n_labels=2):
    """
    Generates fake organ volume data fo use in testing the permutation stats pipeline

    The returned data frame has the followng columns
    index: specimen id
    a series of columns starting with x then label number (eg. x3)
    staging: the staging metric
    line: the line id

    This test tests that the ouput is correct.

    """
    
    wt_organ_vol = pd.read_csv(wt_registration_dir / 'output' / ORGAN_VOLUME_CSV_FILE, index_col=0)
    wt_staging = pd.read_csv(wt_registration_dir / 'output' / STAGING_INFO_FILENAME, index_col=0)
    mut_organ_vol = pd.read_csv(mut_registration_dir / 'output' / ORGAN_VOLUME_CSV_FILE, index_col=0)
    mut_staging = pd.read_csv(mut_registration_dir / 'output' / STAGING_INFO_FILENAME, index_col=0)
    label_meta = target_dir / 'label_info.csv'

    data = run_permutation_stats.prepare_data(wt_organ_vol, wt_staging, mut_organ_vol, mut_staging, label_meta)

    # Check that the staging and line columns are present
    extra_cols = ['staging', 'line']

    flagged_labels = ['2']  # This should not be present in return data as it's flagged

    # Organ labels should be prepended with 'x'. In case we use statsmodels.
    correct_organ_vol_col_names = [f'x{x}' for x in wt_organ_vol if x not in flagged_labels]

    correct_organ_vol_col_names.extend(extra_cols)

    # Check that labels that are flagged to not analyse are removed

    a = set(data.columns)
    b = set(correct_organ_vol_col_names)

    ok_(a == b)



def test_permutation_stats():
    """
    Run the whole permutation based stats pipeline.
    Currently this just checks to see if the pipeline completes without any errors.
    """
    num_perms = 5  # Would do 1000 or more normally
    label_info = registration_root / 'target' / 'label_info.csv'
    label_map = registration_root / 'target' / 'labels.nrrd'
    perm_out_dir = outdir / 'organ_vols_permutation' # Intermediate results go here. Permutation distributions etc.

    run_permutation_stats.run(wt_registration_dir, mut_registration_dir, perm_out_dir, num_perms, log_dependent=False,
                              label_info=label_info, label_map_path=label_map)


@pytest.mark.notest
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


@pytest.mark.notest
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

