"""lama_stats
The main script fo r the non-permutation-based statitics pipelne
This is currently used for the voxel-based data (intensity and jacobians) where permutation testing would be too CPU-intensitve
"""

from pathlib import Path
from logzero import logger as logging

from lama.stats.standard_stats.stats_objects import StatsData
from lama.stats.standard_stats.data_loaders import DataLoader, load_mask
from lama.stats.standard_stats import read_config


def pipeline(data_obj: StatsData):
    """
    Run the stats process on the stast data object

    Parameters
    ----------
    data_obj
        Contains all the input data and we will put the output data there too
    """


def run(config_path: Path, wt_dir: Path, mut_dir: Path, out_dir: Path, target_dir: Path):
    """
    The main function that defines the stats pipeline
    """

    config = read_config.read(config_path)

    mask = load_mask(target_dir, config.mask)
    label_info_file = target_dir / config.get('label_info')

    # Run each data class through the pipeline
    for stats_type in config.stats_types:

        # load the required stats object and data loader
        stats_clas = StatsData.factory(stats_type)
        stats_obj = stats_clas(wt_dir, mut_dir, out_dir, mask, config)
        stats_obj.loader = DataLoader.factory(stats_type)
        stats_obj.load()
        pipeline(stats_obj)





