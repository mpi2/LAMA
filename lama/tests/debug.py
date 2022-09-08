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
from lama.registration_pipeline import run_lama
from lama.tests import (registration_root, mut_registration_dir, wt_registration_dir, test_data_root)


# def delete_previous_files():
#     """
#     Remove the output generated from previous tests. This does not occur directly after the test as we may want to
#     look at the results.
#     """
#     def delete(root: Path):
#         shutil.rmtree(root / 'output', ignore_errors=True)
#         for p in root.iterdir():
#             if str(p).endswith(('.log', 'jobs.csv', 'csv.lock', '.yaml')):
#                 p.unlink()
#
#     delete(wt_registration_dir)
#     delete(mut_registration_dir)


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

    config = test_data_root / 'debugging' /  'debug.toml'

    run_lama.run(config)

    assert True