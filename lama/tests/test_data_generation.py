"""
Test the steps needed to generate wild type and mutant data for use in the statistical analysis

Usage:  pytest -v -m "not notest" test_data_generation.py
The use of -m "not notest" is to be able to omit certain tests with the @pytest.mark.notest decorator
"""

from pathlib import Path
from lama.registration_pipeline import run_lama

import os
import shutil
import pytest

from scripts import lama_job_runner
from . import (registration_root, mut_registration_dir, wt_registration_dir)


@pytest.fixture
def delete_previous_files():
    """
    Remove the output generated from previous tests. This does not occur directly after the test as we may want to
    look at the results.
    """
    def delete(root: Path):
        shutil.rmtree(root / 'output', ignore_errors=True)
        for p in root.iterdir():
            if str(p).endswith(('.log', 'jobs.csv', 'csv.lock', '.yaml')):
                p.unlink()

    delete(wt_registration_dir)
    delete(mut_registration_dir)


def test_make_jobs_file(delete_previous_files):


    config_file = registration_root / 'registration_config.toml'

    lama_job_runner.lama_job_runner(config_file, wt_registration_dir, make_job_file=True)
    lama_job_runner.lama_job_runner(config_file, mut_registration_dir, make_job_file=True)


def test_lama_job_runner():
    """
    Test the lama job runner which was made to utilise multiple machines or the grid.
    This test just uses one machine for the tests at the moment.
    test_make_jobs_file() should run before this to create a jobs file that can be consumed.
    This test should be run before the stats test as it creates data that the stats test needs.

    """

    config_file = registration_root / 'registration_config.toml'

    assert lama_job_runner.lama_job_runner(config_file, wt_registration_dir) is True
    assert lama_job_runner.lama_job_runner(config_file, mut_registration_dir) is True
