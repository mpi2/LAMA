"""lama_stats
The main script fo r the non-permutation-based statitics pipelne
This is currently used for the voxel-based data (intensity and jacobians) where permutation testing would be too CPU-intensitve
"""

from pathlib import Path
from logzero import logger as logging

from lama.stats.standard_stats.stats_objects import Stats
from lama.stats.standard_stats.data_loaders import DataLoader, load_mask
from lama.stats.standard_stats import read_config
from lama.stats.standard_stats import linear_model
from lama.stats.standard_stats import results_writer


def run(config_path: Path,
        wt_dir: Path,
        mut_dir: Path,
        out_dir: Path,
        target_dir: Path):
    """
    The entry point to the stats pipeline.
    Read in the config, and iterate over the stats analysis methods and the mutant lines

    Parameters
    ----------
    config_path
        The toml stats config
    wt_dir
        Root of the wild type data. Should contain mutant line subfolders
    mut_dir
        Root of the mutant data. Should contain mutant line subfolders
    """

    config = read_config.read(config_path)

    mask = load_mask(target_dir, config['mask'])
    label_info_file = target_dir / config.get('label_info')

    # Run each data class through the pipeline
    for stats_type in config['stats_types']:

        # load the required stats object and data loader
        loader_class = DataLoader.factory(stats_type)
        loader = loader_class(wt_dir, mut_dir, mask, config)

        # NOTE: This is where we could parallelise
        for line_input_data in loader.line_iterator():

            stats_class = Stats.factory(stats_type)
            stats_obj = stats_class(line_input_data, stats_type)

            stats_obj.stats_runner = linear_model.lm_r
            stats_obj.run_stats()

            writer = results_writer.factory(stats_type)
            writer(stats_obj, mask, out_dir, stats_type, label_info_file)
