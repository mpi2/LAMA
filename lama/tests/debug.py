"""
For development and debugging, this test can be run. It will use the config debug.toml

Usage:  pytest -qq -m "not notest" debug.py
"""

from pathlib import Path
from lama.registration_pipeline import run_lama
import logzero

import os
import shutil
import pytest
import logging
from lama.scripts import lama_job_runner
from lama.tests import (registration_root, mut_registration_dir, wt_registration_dir)


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


# @pytest.mark.notest
def test_lama_job_runner():
    """
    Test the lama job runner which was made to utilise multiple machines or the grid.
    This test just uses one machine for the tests at the moment.
    test_make_jobs_file() should run before this to create a jobs file that can be consumed.
    This test should be run before the stats test as it creates data that the stats test needs.


    NOTE this test should be at bottom of file as it should be ru last
    The oututs of these tests are consumed by the stats test.
    """

    config = registration_root / 'debug.toml'

    delete_previous_files()

    lama_job_runner.lama_job_runner(config, wt_registration_dir, make_job_file=True, log_level=logging.ERROR)
    lama_job_runner.lama_job_runner(config, wt_registration_dir, log_level=logging.ERROR)

    assert True