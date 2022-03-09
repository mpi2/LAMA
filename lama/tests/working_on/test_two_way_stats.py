from pathlib import Path
from lama.registration_pipeline import run_lama
from lama.scripts import lama_stats
from lama.scripts import lama_job_runner
import logging
import pytest

# Import paths from __init__.py
from lama.tests import (stats_config_dir)

wt_dir = Path(
    "E:/Bl6_data/211014_g_by_back/g_by_back_data/baseline")
mut_dir = Path(
    "E:/Bl6_data/211014_g_by_back/g_by_back_data/mutants")
treat_dir = Path(
    "E:/Bl6_data/211014_g_by_back/g_by_back_data/treatment")
inter_dir = Path(
    "E:/Bl6_data/211014_g_by_back/g_by_back_data/mut_treat")

cfg = Path(
    "E:/Bl6_data/211014_g_by_back/g_by_back_data/generate_data.toml")

stats_cfg = Path(
    "E:/Bl6_data/211014_g_by_back/stats_with_BH_correction/stats.toml")


def test_lama_job_runner():
    """
    Test the lama job runner which was made to utilise multiple machines or the grid.
    This test just uses one machine for the tests at the moment.
    test_make_jobs_file() should run before this to create a jobs file that can be consumed.
    This test should be run before the stats test as it creates data that the stats test needs.
    NOTE this test should be at bottom of file as it should be ru last
    The oututs of these tests are consumed by the stats test.
    """

    print(f"\n{'#' * 8} Doing config {cfg.name} {'#' * 8}")

    lama_job_runner.lama_job_runner(cfg, wt_dir, make_job_file=True, log_level=logging.ERROR)
    lama_job_runner.lama_job_runner(cfg, wt_dir, log_level=logging.ERROR)

    lama_job_runner.lama_job_runner(cfg, mut_dir, make_job_file=True, log_level=logging.ERROR)
    lama_job_runner.lama_job_runner(cfg, mut_dir, log_level=logging.ERROR)

    lama_job_runner.lama_job_runner(cfg, treat_dir, make_job_file=True, log_level=logging.ERROR)
    lama_job_runner.lama_job_runner(cfg, treat_dir, log_level=logging.ERROR)

    lama_job_runner.lama_job_runner(cfg, inter_dir, make_job_file=True, log_level=logging.ERROR)
    lama_job_runner.lama_job_runner(cfg, inter_dir, log_level=logging.ERROR)


def test_g_by_e_reg():
    """
    lama has ony one arg, the config file. Loop over all the configs to test and
    run with lama.
    """

    run_lama.run(cfg)


def test_two_way_stats():
    """
    tests the two_ways_stats component
    Returns
    For each folder:
             intensity (genotype, environment and interaction)
             jacobians (genotype, env etc. )
             organ_volumes (geno, env etc.)
    -------

    """
    lama_stats.run(stats_cfg, wt_dir, mut_dir, treat_dir, inter_dir)
