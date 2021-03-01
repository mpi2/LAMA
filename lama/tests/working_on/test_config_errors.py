"""
Test the steps needed to generate wild type and mutant data for use in the statistical analysis

Usage:  pytest -v -m "not notest" test_data_generation.py
The use of -m "not notest" is to be able to omit certain tests with the @pytest.mark.notest decorator
"""

from pathlib import Path
from lama.registration_pipeline import run_lama, validate_config
from lama.common import cfg_load

import os
import shutil
import pytest

from lama.scripts import lama_job_runner
from lama.tests import (registration_root, mut_registration_dir, wt_registration_dir)

@pytest.mark.notest
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



def test_config_errors():
    """
    Read in the current config that shuld work

    """

    # config_file = registration_root / 'registration_config.toml'
    config_file = registration_root / 'registration_config.toml'
    config = cfg_load(config_file)

    # Staging = embryo_volume needs at least one similarity/affine stage to work
    # for i, stage in enumerate(config['registration_stage_params']):
    #     if stage['elastix_parameters']['Transform'] in ['EulerTransform', 'AffineTransform']:
    #         del(config['registration_stage_params'][i])

    config['registration_stage_params'][:] = [x for x in config['registration_stage_params'] if
                                              x['elastix_parameters']['Transform'] not in ['SimilarityTransform', 'AffineTransform']]

    cfg = validate_config.LamaConfig(config, config_file)

    #
    # assert lama_job_runner.lama_job_runner(config_file, wt_registration_dir) is True
    # assert lama_job_runner.lama_job_runner(config_file, mut_registration_dir) is True
