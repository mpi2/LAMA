import yaml
import logging
import sys
from os.path import join, dirname
import difflib
import enum

sys.path.insert(0, join(dirname(__file__), '..'))
from lib.addict import Dict

TOP_LEVEL_KNOWN_OPTIONS = (
    'data', 'fixed_mask', 'label_map', 'littermate_controls', 'mut_staging_file', 'wt_staging_file', 'n1',
    'voxel_size', 'invert_config_file', 'project_name', 'label_names', 'blur_fwhm', 'formulas', 'use_auto_staging',
    'littermate_pattern', 'mutant_ids', 'inverted_masks'
)

# TOP_LEVEL_KNOWN_OPTIONS = (  # make it an enum
#     'data', 'fixed_mask', 'label_map', 'littermate_controls', 'mut_staging_file', 'wt_staging_file', 'n1',
#     'voxel_size', 'invert_config_file', 'project_name', 'label_names', 'blur_fwhm', 'formulas', 'use_auto_staging',
#     'littermate_pattern', 'mutant_ids'
# )



STATS_JOB_KNOWN_OPTIONS = (
    'wt', 'mut', 'normalisation_roi'
)


def validate(config_path):
    """
    Get the config and check for paths
    """
    try:
        with open(config_path) as fh:
            try:
                config = yaml.load(fh)
            except Exception as e:  # Couldn't catch scanner error from Yaml
                logging.error('Error reading stats yaml file\n\n{}'.format(e))
                sys.exit()
    except IOError:
        logging.error("cannot find or open stats config file: {}".format(config_path))
        sys.exit(1)
    addict_config = Dict(config)

    incorrect = unkown_options(config, TOP_LEVEL_KNOWN_OPTIONS)
    if incorrect:
        sys.exit('incorrect option "{}" in stats yaml file.\nDo you mean {} ?'.format(*incorrect))

    try:
        config['data']
    except KeyError:
        logging.error("stats config file needs a 'data' entry. Are you usin gthe correct config file?")

    else:
        # Check each data entry for correct options (not done yet
        for key, stats_job_config in config['data'].items():
            incorrect = unkown_options(stats_job_config, STATS_JOB_KNOWN_OPTIONS)
            if incorrect:
                sys.exit('incorrect option "{}" in stats yaml file.\nAvaible options are \n\n{}'.format(incorrect[0],
                                                                                                        '\n'.join(
                                                                                                            STATS_JOB_KNOWN_OPTIONS)))
    return addict_config


def unkown_options(config, available_options):
    for param in config:
        if param not in available_options:
            closest_match = difflib.get_close_matches(param, TOP_LEVEL_KNOWN_OPTIONS)
            return param, closest_match
    return False
