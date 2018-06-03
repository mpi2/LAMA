from os.path import join, realpath, dirname, abspath, splitext
import sys
sys.path.insert(0, abspath(join(dirname(__file__), '../..')))

a = sys.path

from nose.tools import assert_equals, nottest

from run_lama import RegistrationPipeline
from elastix.invert import batch_invert_transform_parameters
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

@nottest
def test_lama():
    """
    lama has ony one arg, the config file. Loop over all the configs to test and
    run with lama.
    """

    config_path = abspath(join(baseline_input_dir, lama_configs[0]))
    RegistrationPipeline(config_path)

# @nottest
def test_invert_transforms():
    """
    Test inverting the elastix transform parameters.
    Needs the 'test_lama()' test to have run previously so that the baseline registration data is present
    """

    reg_outdir =  abspath(join(baseline_input_dir, 'output'))
    config_path = abspath(join(baseline_input_dir, lama_configs[0]))
    outdir = join(reg_outdir,'inverted_transforms')
    invert_config = join(outdir, 'invert.yaml')
    batch_invert_transform_parameters(config_path, invert_config, outdir, 1, log=True, noclobber=True)

@nottest
def test_invert_mask():
    pass


@nottest
def test_phenodetect():
    """
    Runs the mutants. Needs the 'test_lama()' test to have run previously so that the baseline
    data is present
    """
    # Config path is one modified form the one that ran the baselines
    phenodetect_config_name = splitext(lama_configs[0])[0] + '_pheno_detect' + '.yaml'
    config_path = abspath(join(baseline_input_dir, phenodetect_config_name))
    PhenoDetect(config_path, mutant_input_dir)