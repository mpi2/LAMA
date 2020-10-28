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

import gc

from lama.common import cfg_load
from lama.stats.standard_stats.stats_objects import Stats, OrganVolume
from lama.stats.standard_stats.data_loaders import DataLoader, load_mask, LineData
from lama.stats.standard_stats.results_writer import ResultsWriter
from lama import common
from lama.stats import linear_model
from lama.elastix.invert_volumes import InvertHeatmap
from lama.img_processing.normalise import Normaliser
from lama.qc import organ_vol_plots


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
        The root output directory. Will be made if not existing

    target_dir
        Contains the population average, masks, label_maps and label infor files
        All Volumes should have been padded to the same size before registration.

    lines_to_process
        list: optional mutant line ids to process only.
        None: process all lines
    """

    if not (wt_dir / 'output').is_dir():
        raise FileNotFoundError(f'{wt_dir / "output"} folder with registration results is not present')
    if not (mut_dir / 'output').is_dir():
        raise FileNotFoundError(f'{mut_dir / "output"} folder with registration results is not present')
    try:
        out_dir.mkdir(exist_ok=True)
    except FileNotFoundError:
        raise FileNotFoundError('Cannot create output folder')

    master_log_file = out_dir / f'{common.date_dhm()}_stats.log'
    logzero.logfile(str(master_log_file))
    logging.info(common.git_log())
    logging.info('### Started stats analysis ###}')

    stats_config = cfg_load(config_path)

    mask = load_mask(target_dir, stats_config['mask'])
    label_info_file = target_dir / stats_config.get('label_info')  # What if not exists
    label_map_file = target_dir / stats_config.get('label_map')
    label_map = common.LoadImage(label_map_file).array

    memmap = stats_config.get('memmap')
    if memmap:
        logging.info('Memory mapping input data')

    baseline_file = stats_config.get('baseline_ids')
    if baseline_file:
        baseline_file = config_path.parent / baseline_file

    mutant_file = stats_config.get('mutant_ids')
    if mutant_file:
        mutant_file = config_path.parent / mutant_file

    # Run each data class through the pipeline.
    for stats_type in stats_config['stats_types']:

        logzero.logfile(str(master_log_file))
        logging.info(f"---Doing {stats_type} analysis---")
        
        gc.collect()
        
        # load the required stats object and data loader
        loader_class = DataLoader.factory(stats_type)

        loader = loader_class(wt_dir, mut_dir, mask, stats_config, label_info_file, lines_to_process=lines_to_process,
                              baseline_file=baseline_file, mutant_file=mutant_file, memmap=memmap)

        # Only affects organ vol loader.
        if not stats_config.get('normalise_organ_vol_to_mask'):
            loader.norm_to_mask_volume_on = False

        # Currently only the intensity stats get normalised
        loader.normaliser = Normaliser.factory(stats_config.get('normalise'), stats_type)  # move this into subclass

        logging.info("Start iterate through lines")
        common.logMemoryUsageInfo()
  
        line_iterator = loader.line_iterator()
        line_input_data = None
 
        while True:
            try:
                line_input_data = next(line_iterator)
                logging.info(f"Data for line {line_input_data.line} loaded")
                common.logMemoryUsageInfo()
         
      
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
      
                logging.info('Statistical analysis finished.')
                common.logMemoryUsageInfo()
                
                logging.info('Writing results...')
                
                rw = ResultsWriter.factory(stats_type)
                writer = rw(stats_obj, mask, line_stats_out_dir, stats_type, label_map, label_info_file)
                
                logging.info('Finished writing results.')
                common.logMemoryUsageInfo()
                #
                # if stats_type == 'organ_volumes':
                #     c_data = {spec: data['t'] for spec, data in stats_obj.specimen_results.items()}
                #     c_df = pd.DataFrame.from_dict(c_data)
                #     # cluster_plots.tsne_on_raw_data(c_df, line_stats_out_dir)
 
      
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
                break;            
         



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
