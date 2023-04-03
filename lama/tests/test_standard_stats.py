"""

Currently just running and making sure there's no uncaught exceptions.
TODO: check that the correct output is generated too

In future we should put this in a web-accessible place

Usage:  pytest test_standard_stats.py
"""

import tempfile
import pytest
import toml
from pathlib import Path

from lama.tests import (wt_registration_dir, mut_registration_dir, target_dir,
                        stats_output_dir)

from lama.stats.standard_stats import lama_stats_new

root_config = dict(
    stats_types=[
        'organ_volumes',
         'intensity',
        'jacobians'
    ],
    normalise_organ_vol_to_mask=True,
    use_log_jacobians = True,
    reg_folder='similarity',
    jac_folder='similarity',
    mask='mask_tight.nrrd',
    label_info='label_info.csv',
    label_map='labels.nrrd',
    blur_fwhm=100,
    voxel_size=14.0,
    invert_stats=True,
    normalise='mask',
    use_staging=True,
    memmap=True
)


@pytest.fixture()
def get_config() -> Path:
    """
    Fixture to get stats config and to write it to temporary file

    Returns
    -------
    location of stats config temporary file
    """
    def make_config(config_updates={}):
        # Create a nested function to allow the fixture to take arguments
        c = root_config.copy()
        c.update(config_updates)

        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.toml') as fh:

            fh.write(toml.dumps(c))

        return c, fh.name

    return make_config





@pytest.mark.skip
def test_all(get_config):
    """
    Run the stats module. The data required for this to work must be initially made
    by running the test:  tests/test_data_generation.py
    """
    config, config_file = get_config()

    lama_stats_new.run(config_file, wt_registration_dir, mut_registration_dir, stats_output_dir, target_dir)

@pytest.mark.skip
def test_no_mask(get_config):
    config, config_file = get_config({'mask': None})

    # Should get value errors if mask not specified
    # with pytest.raises(ValueError):
    #     lama_stats_new.run(config_file, wt_registration_dir, mut_registration_dir, stats_output_dir, target_dir)


# # @nottest
# # def test_no_use():
# #     """
# #     For the organ volume test, if the meta data dataframe has a column of 'no_analysis', do not include the label in the
# #     analysis
# #     """
# #     pass
#
# @nottest
# def test_mutant_id_subset():
#     """
#     Thest whether specifying mutant id subset works
#     """
#
#     config = stats_config_dir / 'new_stats_config_test_mutant_id_subset.toml'
#     lama_stats_new.run(config, wt_registration_dir, mut_registration_dir, stats_output_dir, target_dir)
#
# @nottest
# def test_specifiy_line():
#     """
#     Run the stats module. The data requirted for this to work must be initially made
#     by running the test:  tests/test_lama.py:test_lama_job_runner()
#     """
#     line = 'mutant_1'  # Only run this line
#
#     # Remove any previous output
#     if stats_output_dir.is_dir():
#         shutil.rmtree(stats_output_dir)
#
#     config = stats_config_dir / 'new_stats_config.toml'
#     lama_stats_new.run(config, wt_registration_dir, mut_registration_dir, stats_output_dir, target_dir, lines_to_process=[line])
#
#     output_lines = [x.name for x in stats_output_dir.iterdir() if x.is_dir()]
#     assert_equal([line],  output_lines)
#
#
# @nottest
# def test_erroneuos_configs():
#     error_config_dir = stats_config_dir / 'erroneous_configs'
#
#     for config in error_config_dir.iterdir():
#         lama_stats_new.run(config, wt_registration_dir, mut_registration_dir, stats_output_dir, target_dir)
#
# # @nottest
# # def test_organ_vols():
# #     config = join(CONFIG_DIR, 'organ_vols.yaml')
# #     run_lama_stats.run(config)
# #
# # @nottest
# # def test_glcm():
# #     config = join(CONFIG_DIR, 'test_glcm.yaml')
# #     run_lama_stats.run(config)
# #
# # @nottest
# # def test_from_arbitrary_directory():
# #     config = join(CONFIG_DIR, 'all_specimens.yaml')
# #     os.chdir(os.path.expanduser('~'))
# #     run_lama_stats.run(config)
# #
# # @nottest
# # def test_missing_config():
# #     config = join(CONFIG_DIR, 'all_specimens.flaml')
# #     assert_raises(Exception, run_lama_stats.run, config)
# #
# # @nottest
# # def test_misc():
# #     config = join(CONFIG_DIR, 'misc_test.yaml')
# #     run_lama_stats.run(config)
# #
# #
# # @nottest
# # def test_incorrect_staging_file():
# #     pass
#
#
