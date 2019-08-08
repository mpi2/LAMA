"""lama_stats
The main script for the non-permutation-based statitics pipelne
This is currently used for the voxel-based data (intensity and jacobians) where permutation testing would be too
CPU-intensive

Outline of the stats pipeline
------------------------------
To do 060819


"""

from pathlib import Path
from typing import Union, List

from logzero import logger as logging
import logzero

from lama.stats.standard_stats.stats_objects import Stats
from lama.stats.standard_stats.data_loaders import DataLoader, load_mask, LineData
from lama.stats.standard_stats import read_config
from lama.stats.standard_stats import linear_model
from lama.stats.standard_stats.results_writer import ResultsWriter
from lama import common
from lama.stats import cluster_plots
from lama.elastix.invert_volumes import InvertHeatmap
from lama.img_processing.normalise import Normaliser


def run(config_path: Path,
        wt_dir: Path,
        mut_dir: Path,
        out_dir: Path,
        target_dir: Path,
        lines_to_process: Union[List, None] = None
        ):
    """
    The entry point to the stats pipeline.
    Read in the stats_config, and iterate over the stats analysis methods and the mutant lines

    Parameters
    ----------
    config_path
        The lama stats_config (in TOML format)

    wt_dir
        Root of the wild type data. Should contain mutant line subfolders

    mut_dir
        Root of the mutant data. Should contain mutant line subfolders

    out_dir
        The root utput directory

    target_dir
        Contains the population average, masks, label_maps and label infor files
        All Volumes should have been padded to the same size before registration.

    lines_to_process
        list: optional mutant line ids to process only.
        None: process all lines
    """
    out_dir.mkdir(exist_ok=True)
    master_log_file = out_dir / f'{common.date_dhm()}_stats.log'
    logzero.logfile(str(master_log_file))
    logging.info('### Started stats analysis ###}')

    stats_config = read_config.read(config_path)

    mask = load_mask(target_dir, stats_config['mask'])
    label_info_file = target_dir / stats_config.get('label_info')  # What if not exists
    label_map_file = target_dir / stats_config.get('label_map')
    label_map = common.LoadImage(label_map_file).array

    baseline_file = stats_config.get('baseline_ids')
    if baseline_file:
        baseline_file = config_path.parent / baseline_file

    mutant_file = stats_config.get('mutant_ids')
    if mutant_file:
        mutant_file = config_path.parent / mutant_file

    # Run each data class through the pipeline.
    for stats_type in stats_config['stats_types']:
        logzero.logfile(str(master_log_file))
        logging.info(f"Doing {stats_type} analysis")
        # load the required stats object and data loader
        loader_class = DataLoader.factory(stats_type)

        loader = loader_class(wt_dir, mut_dir, mask, stats_config, label_info_file, lines_to_process=lines_to_process,
                              baseline_file=baseline_file, mutant_file=mutant_file)

        loader.normaliser = Normaliser.factory(stats_config.get('normalise'), stats_type)  # move this into subclass

        for line_input_data in loader.line_iterator():  # NOTE: This might be where we could parallelise

            line_id = line_input_data.line

            line_stats_out_dir = out_dir / line_id / stats_type

            line_stats_out_dir.mkdir(parents=True, exist_ok=True)
            line_log_file = line_stats_out_dir / f'{common.date_dhm()}_stats.log'
            logzero.logfile(str(line_log_file))

            logging.info(f"Processing line: {line_id}")

            stats_class = Stats.factory(stats_type)
            stats_obj = stats_class(line_input_data, stats_type, stats_config.get('use_staging', True))

            stats_obj.stats_runner = linear_model.lm_r
            stats_obj.run_stats()

            logging.info('statistical analysis finished. Writing results. ')

            rw = ResultsWriter.factory(stats_type)
            writer = rw(stats_obj, mask, line_stats_out_dir, stats_type, label_map, label_info_file)

            # cluster_plots.tsne_on_raw_data(stats_obj.cluster_data(), line_stats_out_dir)

            if stats_config.get('invert_stats'):
                if writer.line_heatmap:  # Organ vols wil not have this
                    # How do I now sensibily get the path to the invert.yaml
                    # get the invert_configs for each specimen in the line
                    line_heatmap = writer.line_heatmap
                    line_reg_dir = mut_dir / 'output' / line_id
                    invert_heatmaps(line_heatmap, line_stats_out_dir, line_reg_dir, line_input_data)


def invert_heatmaps(heatmap: Path,
                    stats_outdir: Path,
                    reg_outdir: Path,
                    input_: LineData):
    """
    Invert the stats heatmaps from a single line back onto inputs or registered volumes

    Parameters
    ----------
    line_dir
        The registration output directory for a line
    input_
        Has paths for data locations
    outdir
        Where to put the inverted heatmaps

    Returns
    -------

    """
    #  Do some logging
    inverted_heatmap_dir = stats_outdir / 'inverted_heatmaps'
    common.mkdir_force(inverted_heatmap_dir)

    for spec_id in input_.mutant_ids():
        # Should not have to specify the path to the inv config again
        invert_config = reg_outdir /  spec_id/ 'output' / 'inverted_transforms' / 'invert.yaml'

        inv = InvertHeatmap(invert_config, heatmap, inverted_heatmap_dir)
        inv.run()
