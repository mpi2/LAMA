
"""
These functions test the lama registration pipeline
"""

from os.path import join, realpath, dirname, abspath, splitext
from pathlib import Path
import sys

current_dir = dirname(realpath(__file__))
sys.path.insert(0, abspath(join(dirname(__file__), '../..')))

import warnings
warnings.filterwarnings('ignore')

from nose.tools import nottest, assert_raises

from run_lama import RegistrationPipeline
from job_runner import lama_job_runner
from elastix.invert import batch_invert_transform_parameters, InvertLabelMap
from pheno_detect import PhenoDetect

test_data_root = join(current_dir, '..', 'test_data', 'registration_test_data')

INPUT_DIR = join(test_data_root, 'input_data')

current_dir = dirname(realpath(__file__))

baseline_input_dir = join(test_data_root, 'baseline')

mutant_input_dir = join(test_data_root, 'mutant')

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


@nottest
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
def test_invert_labels():
    """
    Test inverting a labelmap using the inverted transform parameters
    Needs test_invert_transforms to have been run previously

    Returns
    -------
    """
    labelmap_path = abspath(join(baseline_input_dir, 'target', 'v23_merge_recut_cleaned_flipped.nrrd'))
    invert_config = abspath(join(baseline_input_dir, 'output', 'inverted_transforms', 'invert.yaml'))
    outdir = abspath(join(baseline_input_dir, 'output', 'inverted_lables'))
    inv = InvertLabelMap(invert_config, labelmap_path, outdir, threads=1, noclobber=False)
    inv.run()


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


# @nottest
def test_lama_job_runner_baselines():
    """
    Testing lama job runner for baselines

    This is a work in progress

    Test the lama job runner which was made to utilise multiple machines or the grid.
    Just using one mahine for the tests at the moment. Add multiple machines later
    -------

    """

    root_folder = Path(baseline_input_dir)

    config_file = Path(test_data_root) / 'registration_config.yaml'

    assert_raises(SystemExit, lama_job_runner, config_file, root_folder, type_='registration')


# @nottest
def test_lama_job_runner_mutants():

    # Run the mutants as well. This will give data we can use for the stats
    root_folder = Path(mutant_input_dir)

    config_file = Path(test_data_root) / 'registration_config.yaml'

    assert_raises(SystemExit, lama_job_runner, config_file, root_folder, type_='registration')

