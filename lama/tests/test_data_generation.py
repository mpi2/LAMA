"""
Test the steps needed to generate wild type and mutant data for use in the statistical analysis

Usage:  pytest -qq -m "not notest" test_data_generation.py
The use of -m "not notest" is to be able to omit certain tests with the @pytest.mark.notest decorator
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

    configs = registration_root.glob('*.toml')

    for cfg in configs:
        delete_previous_files()

        print(f"\n{'#'*8} Doing config {cfg.name} {'#'*8}")

        lama_job_runner.lama_job_runner(cfg, wt_registration_dir, make_job_file=True, log_level=logging.ERROR)
        lama_job_runner.lama_job_runner(cfg, wt_registration_dir, log_level=logging.ERROR)

        lama_job_runner.lama_job_runner(cfg, mut_registration_dir, make_job_file=True, log_level=logging.ERROR)
        lama_job_runner.lama_job_runner(cfg, mut_registration_dir, log_level=logging.ERROR)



#TODO
# QC subset test



#
# @pytest.mark.notest
# def test_make_jobs_file():
#
#
#     config_file = registration_root / 'registration_config.toml'
#
#     lama_job_runner.lama_job_runner(config_file, wt_registration_dir, make_job_file=True)
#     lama_job_runner.lama_job_runner(config_file, mut_registration_dir, make_job_file=True)
#
# @pytest.mark.notest
# def test_lama_job_runner_reverse_reg_only():
#     """
#     Tests out doing only the propagation of atlas to input images. The first stage of the noraml foraward registration
#     of input the pop abg is carried out. This creates the rigid registered input which is then used as target to
#     map atlas to.
#     """
#     config_file = registration_root / 'registration_config_reverse_reg_only.toml'
#     assert lama_job_runner.lama_job_runner(config_file, mut_registration_dir) is True
#
# @pytest.mark.notest
# def test_lama_job_runner_pyramid():
#     """
#     map atlas to.
#     """
#     config_file = registration_root / 'registration_config.toml'
#     assert lama_job_runner.lama_job_runner(config_file, mut_registration_dir) is True
#
#
# @pytest.mark.notest
# def test_lama_reg():
#     """
#     Test using the lama registration script without jaob runner wrapepr
#     Returns
#     -------
#
#     """
#     # delete_previous_files()
#     config_file = registration_root / 'registration_config.toml'
#     # Needs to be in same folder as inputs
#     dest = wt_registration_dir / 'registration_config.toml'
#     # shutil.copyfile(config_file, dest)
#     assert run_lama.run(dest) is True
#
#
#
#
# @pytest.mark.notest
# def test_lama_job_runner_secondary_segmentation():
#     """
#     Tests out doing only the propagation of atlas to input images. The first stage of the noraml foraward registration
#     of input the pop abg is carried out. This creates the rigid registered input which is then used as target to
#     map atlas to.
#     """
#     config_file = registration_root / 'registration_config_reverse_reg_only.toml'
#     assert lama_job_runner.lama_job_runner(config_file, mut_registration_dir) is True


