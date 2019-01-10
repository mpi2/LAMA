"""lama_stats
The main script for the non-permutation-based statitics pipelne
This is currently used for the voxel-based data (intensity and jacobians) where permutation testing would be too CPU-intensive \
or wher
"""

from pathlib import Path
from logzero import logger as logging
import logzero

from lama.stats.standard_stats.stats_objects import Stats
from lama.stats.standard_stats.data_loaders import DataLoader, load_mask, LineData
from lama.stats.standard_stats import read_config
from lama.stats.standard_stats import linear_model
from lama.stats.standard_stats.results_writer import ResultsWriter
from lama import common
from lama.stats import cluster_plots
from lama.elastix.invert_volumes import InvertStats
from lama.registration_pipeline.validate_config import LamaConfig


def run(config_path: Path,
        wt_dir: Path,
        mut_dir: Path,
        out_dir: Path,
        target_dir: Path,
        # This is just a test. Maybe passing in the orginal lama config is the way to go so we don't have to work out all the paths agaain
        lama_config: Path
        ):
    """
    The entry point to the stats pipeline.
    Read in the stats_config, and iterate over the stats analysis methods and the mutant lines

    Parameters
    ----------
    config_path
        The lama stats_config (in TOML format). This stats_config should be in the orginal folder it was used to generate the
        registration data.
    wt_dir
        Root of the wild type data. Should contain mutant line subfolders
    mut_dir
        Root of the mutant data. Should contain mutant line subfolders
    out_dir
        The root utput directory
    target_dir
        Contains the population average, masks, label_maps and label infor files
        All Volumes should have been padded to the same size before registration.
    """

    config = LamaConfig(lama_config)

    master_log_file = out_dir / f'{common.date_dhm()}_stats.log'
    logzero.logfile(master_log_file)
    logging.info('### Started stats analysis ###}')

    stats_config = read_config.read(config_path)

    mask = load_mask(target_dir, stats_config['mask'])
    label_info_file = target_dir / stats_config.get('label_info')

    # Run each data class through the pipeline.
    for stats_type in stats_config['stats_types']:

        # load the required stats object and data loader
        loader_class = DataLoader.factory(stats_type)
        loader = loader_class(wt_dir, mut_dir, mask, stats_config)

        for line_input_data in loader.line_iterator():  # NOTE: This might be where we could parallelise

            line_stats_out_dir = out_dir / line_input_data.line / stats_type
            line_stats_out_dir.mkdir(parents=True, exist_ok=True)

            stats_class = Stats.factory(stats_type)
            stats_obj = stats_class(line_input_data, stats_type)

            stats_obj.stats_runner = linear_model.lm_r
            stats_obj.run_stats()

            writer = ResultsWriter.factory(stats_type)
            writer(stats_obj, mask, line_stats_out_dir, stats_type, label_info_file)

            cluster_plots.tsne_on_raw_data(line_input_data, line_stats_out_dir)

            if stats_config.get('invert_stats'):
                # How do I now sensibily get the path to the invert.yaml
                # get the invert_configs for each specimen in the line
                inv_configs = get_inv_configs()
                inv = InvertStats(config_path, invertable, outdir)
                inv.run()


            # results_writer.pvalue_fdr_plot(stats_obj, )

def get_inv_configs():
    pass