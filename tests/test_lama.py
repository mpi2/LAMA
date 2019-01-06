
"""
These functions test the lama registration pipeline

Usage:

"""

from os.path import join, realpath, dirname

from pathlib import Path

current_dir = dirname(realpath(__file__))

import warnings
warnings.filterwarnings('ignore')

from nose.tools import nottest, assert_raises

from lama.registration_pipeline.run_lama import run
from lama.registration_pipeline.validate_config import LamaConfig
from lama.elastix import deformations
from lama.elastix import invert_volumes
from lama.registration_pipeline.run_lama import invert_volumes
from lama.elastix.invert_transforms import batch_invert_transform_parameters


from scripts import lama_job_runner  # yes import

test_data_root = join(current_dir, 'test_data')

registration_root = join(test_data_root, 'registration_test_data')

INPUT_DIR = join(test_data_root, 'input_data')

current_dir = dirname(realpath(__file__))

baseline_input_dir = join(registration_root, 'baseline')
population_test_dir = join(test_data_root, 'population_average_data')



mutant_input_dir = join(test_data_root, 'mutant')

lama_configs = [
    'lama.yaml'
]


@nottest
def test_lama_job_runner():
    """
    Testing lama job runner for baselines

    This is a work in progress

    Test the lama job runner which was made to utilise multiple machines or the grid.
    Just using one mahine for the tests at the moment. Add multiple machines later
    -------

    """

    root_folder = Path(registration_root)/ 'baseline'

    config_file = Path(registration_root) / 'registration_config.toml'

    assert_raises(SystemExit, lama_job_runner.lama_job_runner, config_file, root_folder)

    root_folder = Path(registration_root)/ 'mutant'

    assert_raises(SystemExit, lama_job_runner.lama_job_runner, config_file, root_folder)


@nottest
def test_population_average():
    """
    lama has ony one arg, the config file. Loop over all the configs to test and
    run with lama.
    """

    config_file = Path(population_test_dir)/ 'registration_config.toml'
    run(config_file)

@nottest
def test_deformations():
    config_file = Path(population_test_dir) / 'registration_config.toml'
    deformations.make_deformations_at_different_scales(config_file)

@nottest
def test_make_inversion_files():
    config_file = Path(population_test_dir) / 'registration_config.toml'
    invert_volumes.batch_invert_transform_parameters(config_file)

def test_invert_transforms():
    config_file = Path(population_test_dir) / 'registration_config.toml'
    config = LamaConfig(config_file)
    batch_invert_transform_parameters(config)

@nottest
def test_invert_volumes():
    config_file = Path(population_test_dir) / 'registration_config.toml'
    config = LamaConfig(config_file)
    invert_volumes(config)

@nottest
def test_validate_config():
    config_file = Path(population_test_dir) / 'registration_config.toml'
    run(config_file)

# @nottest
# def test_invert_transforms():
#     """
#     Test inverting the elastix transform parameters.
#     Needs the 'test_lama()' test to have run previously so that the baseline registration data is present
#     """
#
#     reg_outdir =  abspath(join(baseline_input_dir, 'output'))
#     config_path = abspath(join(baseline_input_dir, lama_configs[0]))
#     outdir = join(reg_outdir,'inverted_transforms')
#     invert_config = join(outdir, 'invert.yaml')
#     batch_invert_transform_parameters(config_path, invert_config, outdir, 1, log=True, noclobber=True)
#
#
# @nottest
# def test_invert_labels():
#     """
#     Test inverting a labelmap using the inverted transform parameters
#     Needs test_invert_transforms to have been run previously
#
#     Returns
#     -------
#     """
#     labelmap_path = abspath(join(baseline_input_dir, 'target', 'v23_merge_recut_cleaned_flipped.nrrd'))
#     invert_config = abspath(join(baseline_input_dir, 'output', 'inverted_transforms', 'invert.yaml'))
#     outdir = abspath(join(baseline_input_dir, 'output', 'inverted_lables'))
#     inv = InvertLabelMap(invert_config, labelmap_path, outdir, threads=1, noclobber=False)
#     inv.run()
#
#
# @nottest
# def test_phenodetect():
#     """
#     Runs the mutants. Needs the 'test_lama()' test to have run previously so that the baseline
#     data is present
#     """
#     # Config path is one modified form the one that ran the baselines
#     phenodetect_config_name = splitext(lama_configs[0])[0] + '_pheno_detect' + '.yaml'
#     config_path = abspath(join(baseline_input_dir, phenodetect_config_name))
#     PhenoDetect(config_path, mutant_input_dir)
#
#
#
#
# def test_job_runner_cli():
#     pass
#
#
# @nottest
# def test_lama_job_runner_mutants():
#
#     # Run the mutants as well. This will give data we can use for the stats
#     root_folder = Path(mutant_input_dir)
#
#     config_file = Path(test_data_root) / 'registration_config.yaml'
#
#     assert_raises(SystemExit, lama_job_runner, config_file, root_folder, type_='registration')



