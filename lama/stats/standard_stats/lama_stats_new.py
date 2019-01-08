"""lama_stats
The main script for the non-permutation-based statitics pipelne
This is currently used for the voxel-based data (intensity and jacobians) where permutation testing would be too CPU-intensive \
or wher
"""

from pathlib import Path
from logzero import logger as logging
import logzero

from lama.stats.standard_stats.stats_objects import Stats
from lama.stats.standard_stats.data_loaders import DataLoader, load_mask, InputData
from lama.stats.standard_stats import read_config
from lama.stats.standard_stats import linear_model
from lama.stats.standard_stats.results_writer import ResultsWriter
from lama import common
from lama.stats import cluster_plots


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
    out_dir
        The root utput directory
    target_dir
        Contains the population average, masks, label_maps and label infor files
        All Volumes should have been padded to the same size before registration.
    """

    master_log_file = out_dir / f'{common.date_dhm()}_stats.log'
    logzero.logfile(master_log_file)
    logging.info('### Started stats analysis ###}')

    config = read_config.read(config_path)

    mask = load_mask(target_dir, config['mask'])
    label_info_file = target_dir / config.get('label_info')

    # Run each data class through the pipeline.
    for stats_type in config['stats_types']:

        # load the required stats object and data loader
        loader_class = DataLoader.factory(stats_type)
        loader = loader_class(wt_dir, mut_dir, mask, config)

        for line_input_data in loader.line_iterator():  # NOTE: This might be where we could parallelise

            stats_class = Stats.factory(stats_type)
            stats_obj = stats_class(line_input_data, stats_type)

            stats_obj.stats_runner = linear_model.lm_r
            stats_obj.run_stats()

            writer = ResultsWriter.factory(stats_type)
            w = writer(stats_obj, mask, out_dir, stats_type, label_info_file)

            # This is a bodge for now until I work out how the cluster plots get written
            out_dir = w.out_dir
            cluster_plots.tsne_on_raw_data(line_input_data, out_dir)



            # results_writer.pvalue_fdr_plot(stats_obj, )

def _log_input_data(in_: InputData, stats_type: str):

    logging.info(f'Started stats analysis\nline:{in_.line}\nstats type: {stats_type}')
    logging.info('Using wild type paths\n: {}'.format("\n".join(in_.paths[0])))
    logging.info('Using mutant paths\n: {}'.format("\n".join(in_.paths[1])))
