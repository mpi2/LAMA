"""lama_stats
The main script fo r the non-permutation-based statitics pipelne
This is currently used for the voxel-based data (intensity and jacobians) where permutation testing would be too CPU-intensitve
"""

from pathlib import Path

from lama.stats.standard_stats.stats_objects import StatsData
from lama.stats.standard_stats import read_config


def run(config_path: Path):
    """
    The main function that defines the stats pipeline
    """

    config = read_config.read(config_path)
    stats_data_objs = get_data_(config.stats_types)

def get_data_():
    # Read in the config not sure what that will look like yet
    data_store = StatsData.factory('intensity')

