
"""
These functions test the lama registration pipeline

Usage:  pytest -v -m "not notest" test_population_average.py
The use of -m "not notest" is to be able to omit certain tests with the @pytest.mark.notest decorator

"""


from pathlib import Path
from lama.registration_pipeline import run_lama

import pytest


# Import paths from __init__.py
from lama.tests import (population_test_dir)

@pytest.mark.skip
def test_population_average():
    """
    lama has ony one arg, the config file. Loop over all the configs to test and
    run with lama.
    """

    config_file = Path(population_test_dir)/ 'population_average_config.toml'
    run_lama.run(config_file)



def test_pairwie_average():
    config_file = Path(population_test_dir)/ 'pairwise_population_average_config.toml'
    run_lama.run(config_file)