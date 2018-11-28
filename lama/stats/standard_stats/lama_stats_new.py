"""lama_stats
The main script fo r the non-permutation-based statitics pipelne
This is currently used for the voxel-based data (intensity and jacobians) where permutation testing would be too CPU-intensitve
"""

from pathlib import Path

from lama.stats.standard_stats.stats_objects import StatsData
from lama.stats.standard_stats.data_loaders import DataLoader
from lama.stats.standard_stats import read_config


def pipeline(data_obj: StatsData):
    pass


def run(config_path: Path, wt_dir: Path, mut_dir: Path):
    """
    The main function that defines the stats pipeline
    """

    config = read_config.read(config_path)
    root_out_dir = config_path.parent()

    # Run each data class through the pipeline
    for stats_type in config.stats_types:

        # load the required stats object and data loader
        stats_obj = StatsData.factory(stats_type, wt_dir, mut_dir, root_out_dir)
        stats_obj.loader = DataLoader.factory(stats_type, config)





