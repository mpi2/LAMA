from pathlib import Path
from lama.registration_pipeline import run_lama
from lama.scripts import lama_stats
from lama.scripts import lama_job_runner
import logging
import pytest
from lama.stats.standard_stats.data_loaders import DataLoader, load_mask, LineData, JacobianDataLoader
from lama import common
from lama.img_processing.normalise import Normaliser
import logzero

from lama.stats.standard_stats.results_writer import ResultsWriter
from lama.stats.standard_stats.stats_objects import Stats, OrganVolume
from lama.stats.standard_stats.lama_stats_new import invert_heatmaps

from lama.stats import linear_model

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

target_dir =  Path(
    "E:/Bl6_data/211014_g_by_back/target")

stats_output =  Path("E:/Bl6_data/211014_g_by_back/stats_with_BH_correction/stats_output")

lines_to_process = None

@pytest.mark.skip
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

@pytest.mark.skip
def test_g_by_e_reg():
    """
    lama has ony one arg, the config file. Loop over all the configs to test and
    run with lama.
    """

    run_lama.run(cfg)


def test_two_way_intensities():
    stats_config = common.cfg_load(stats_cfg)

    loader_class = DataLoader.factory("intensity")


    mask = load_mask(target_dir, stats_config['mask'])
    label_info_file = target_dir / stats_config.get('label_info')  # What if not exists
    label_map_file = target_dir / stats_config.get('label_map')
    label_map = common.LoadImage(label_map_file).array

    memmap = stats_config.get('memmap')
    if memmap:
        logging.info('Memory mapping input data')

    baseline_file = stats_config.get('baseline_ids')
    if baseline_file:
        baseline_file = stats_cfg.parent / baseline_file

    mutant_file = stats_config.get('mutant_ids')
    if mutant_file:
        mutant_file = stats_cfg.parent / mutant_file

    loader = loader_class(wt_dir, mut_dir, mask, stats_config, label_info_file, lines_to_process=lines_to_process,
                          baseline_file=baseline_file, mutant_file=mutant_file, memmap=memmap,
                          treatment_dir=treat_dir, interaction_dir=inter_dir)

    loader.normaliser = Normaliser.factory(stats_config.get('normalise'), "intensity") # move this into subclass

    line_iterator = loader.two_way_iterator()
    line_input_data = None

    while True:
        try:
            line_input_data = next(line_iterator)
            logging.info(f"Data for line {line_input_data.line} loaded")
            common.logMemoryUsageInfo()

            line_id = line_input_data.line

            line_stats_out_dir = stats_output / line_id / "intensity"

            line_stats_out_dir.mkdir(parents=True, exist_ok=True)
            line_log_file = line_stats_out_dir / f'{common.date_dhm()}_stats.log'
            logzero.logfile(str(line_log_file))

            stats_class = Stats.factory("intensity")
            stats_obj = stats_class(line_input_data, "intensity", stats_config.get('use_staging', True),
                                    stats_config.get('two_way', False))

            stats_obj.stats_runner = linear_model.lm_r
            stats_obj.run_stats()

            logging.info('Statistical analysis finished.')
            common.logMemoryUsageInfo()

            logging.info('Writing results...')

            rw = ResultsWriter.factory("intensity")
            writer = rw(stats_obj, mask, line_stats_out_dir, "intensity", label_map, label_info_file,
                        stats_config.get('two_way', False))

            if stats_config.get('invert_stats'):
                if writer.line_heatmap:  # Organ vols wil not have this
                    # How do I now sensibily get the path to the invert.yaml
                    # get the invert_configs for each specimen in the line
                    logging.info('Writing heatmaps...')
                    logging.info('Propogating the heatmaps back onto the input images ')
                    line_heatmap = writer.line_heatmap
                    line_reg_dir = mut_dir / 'output' / line_id
                    invert_heatmaps(line_heatmap, line_stats_out_dir, line_reg_dir, line_input_data)
                    logging.info('Finished writing heatmaps.')

            logging.info(f"Finished processing line: {line_id} - All done")
            common.logMemoryUsageInfo()

        except StopIteration:
            if (line_input_data != None):
                logging.info(f"Finish iterate through lines")
                line_input_data.cleanup()
                common.logMemoryUsageInfo()

            break






@pytest.mark.skip
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
    lama_stats.run(stats_cfg, wt_dir, mut_dir, stats_output, treat_dir, inter_dir)
