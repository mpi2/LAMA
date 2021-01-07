"""
Test the permutation-based stats pipeline

Usage
-----

These functions test the lama registration pipeline permutation stats module

Usage:  pytest -q -x -s  -m "not notest"  --tb=short test_run_permutation_stats.py
"""

import shutil

import pytest

from lama.stats.permutation_stats import run_permutation_stats
from lama.tests import (test_data_root, registration_root, wt_registration_dir, mut_registration_dir,
                        permutation_stats_dir, qc_flags_dir)
from lama.common import LamaDataException


@pytest.fixture(scope="session", autouse=True)
def remove_previous_output():
    """
    Remove previous output from this test module. We do not do this after each test as manual inspection of the output
    may be necessary. this runs before all other tests.
    """
    shutil.rmtree(permutation_stats_dir, ignore_errors=True)
    permutation_stats_dir.mkdir()


# @pytest.mark.notest
def test_permutation_stats():
    """
    Run the whole permutation based stats pipeline.
    Currently this just checks to see if the pipeline completes without any errors.
    """
    outdir = permutation_stats_dir / 'output'
    outdir.mkdir()
    num_perms = 5  # Would do 1000 or more normally
    label_info = registration_root / 'target' / 'label_info.csv'
    label_map = registration_root / 'target' / 'labels.nrrd'

    run_permutation_stats.run(wt_registration_dir, mut_registration_dir, outdir, num_perms,
                              label_info=label_info, label_map_path=label_map)


def test_permutation_stats_with_qc_flaggs():
    """
    Run the permutations stats but include a specimen/organ-level qc file to exclude qc-flagged organs
    """
    num_perms = 5  # Would do 1000 or more normally
    label_info = registration_root / 'target' / 'label_info.csv'
    label_map = registration_root / 'target' / 'labels.nrrd'

    for qc_file in qc_flags_dir.iterdir():

        out_dir = permutation_stats_dir / qc_file.stem  # Intermediate results go here. Permutation distributions etc.
        out_dir.mkdir()

        if 'error' in str(qc_file):
            # These qc flag files with errors in should raise a LamaDataError
            with pytest.raises(LamaDataException):
                run_permutation_stats.run(wt_registration_dir, mut_registration_dir, out_dir, num_perms,
                                          label_info=label_info, label_map_path=label_map, qc_file=qc_file)

        else:
            run_permutation_stats.run(wt_registration_dir, mut_registration_dir, out_dir, num_perms,
                                  label_info=label_info, label_map_path=label_map, qc_file=qc_file)


# @pytest.mark.notest
# def test_p_thresholds():
#     """
#     Testing the p_thresholds calculation
#     -------
#     TODO: Add more tests for different cases
#     """
#     import copy
#     import pandas as pd
#
#     # These simulate null distributions from 100 baselines for two organs
#     null = pd.DataFrame.from_dict({
#         '1': np.linspace(0.01, 1, 100),
#         '2': np.linspace(0.01, 1, 100)})
#
#     # These simulate alternative distributions for two organs from 40 lines
#     alt = pd.DataFrame.from_dict({
#         '1': np.linspace(0.00000001, 0.1, 40),
#         '2': np.linspace(0.9, 1, 40)})
#
#     thresh = p_thresholds.get_thresholds(null, alt)
#
#     assert thresh.loc[1, 'p_thresh'] == 0.02  # Gives a p-value threshold of 0.02
#     assert thresh.loc[2, 'p_thresh'] == 1.0    # Gives a p-value threshold of 1.0 as there are no low p-values in the alt distribution


# @pytest.mark.notest
# def test_annotate():
#     # Lines
#     alt_file = Path('/home/neil/git/lama/tests/test_data/stats_test_data/test_output/organ_vols_permutation/alt_line_dist_pvalues.csv')
#     thresholds_file = Path('/home/neil/git/lama/tests/test_data/stats_test_data/test_output/organ_vols_permutation/line_organ_p_thresholds.csv')
#     mutant_dir = Path('/home/neil/git/lama/tests/test_data/registration_test_data/mutant')
#
#     thresholds = pd.read_csv(thresholds_file, index_col=0)
#     alt = pd.read_csv(alt_file, index_col=0)
#
#     run_permutation_stats.annotate(thresholds, alt, mutant_dir)
#
#     # # Specimens
#     alt_file = Path('/home/neil/git/lama/tests/test_data/stats_test_data/test_output/organ_vols_permutation/alt_specimen_dist_pvalues.csv')
#     thresholds_file = Path('/home/neil/git/lama/tests/test_data/stats_test_data/test_output/organ_vols_permutation/specimen_organ_p_thresholds.csv')
#     mutant_dir = Path('/home/neil/git/lama/tests/test_data/registration_test_data/mutant')
#
#     thresholds = pd.read_csv(thresholds_file, index_col=0)
#     alt = pd.read_csv(alt_file, index_col=0)
#
#     run_permutation_stats.annotate(thresholds, alt, mutant_dir)

