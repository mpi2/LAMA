from os.path import join, realpath, dirname, abspath, splitext
import sys
sys.path.insert(0, join(dirname(__file__), '..'))

from nose.tools import assert_equals, nottest

from lama import RegistraionPipeline
from pheno_detect import PhenoDetect
from . import INPUT_DIR

"""
These functions test the lama registration pipeline
"""

current_dir = dirname(realpath(__file__))
baseline_input_dir = join(INPUT_DIR, 'baselines')
mutant_input_dir = join(INPUT_DIR, 'mutant')

lama_configs = [
    'lama.yaml'
]

# @nottest
def test_all_lama_configs():
    """
    lama has ony one arg, the config file. Loop over all the configs to test and
    run with lama.
    """

    config_path = abspath(join(baseline_input_dir, lama_configs[0]))
    RegistraionPipeline(config_path)

@nottest
def test_phenodetect():
    """
    Runs the mutants. Needs the lama test to have run previously so that the baseline
    data is present
    """
    # Config path is one modified form the one that ran the baselines
    phenodetect_config_name = splitext(lama_configs[0])[0] + '_pheno_detect' + '.yaml'
    config_path = abspath(join(baseline_input_dir, phenodetect_config_name))
    PhenoDetect(config_path, mutant_input_dir)