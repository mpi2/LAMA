import yaml
import logging
import sys
from os.path import join, dirname
import os
import difflib
import enum

sys.path.insert(0, join(dirname(__file__), '..'))
from lib.addict import Dict

TOP_LEVEL_KNOWN_OPTIONS = (
    'data', 'littermate_controls', 'n1',
    'voxel_size', 'project_name', 'blur_fwhm', 'formulas', 'use_auto_staging',
    'littermate_pattern', 'mutant_ids', 'wildtype_ids'
)

TOP_LEVEL_KNOWN_PATHS = (
    'output_dir', 'fixed_mask', 'label_map', 'mut_staging_file', 'wt_staging_file',
    'invert_config_file', 'label_names', 'inverted_masks'
)

STATS_JOB_KNOWN_OPTIONS = (
    'wt', 'mut', 'normalisation', 'wt_inverted_masks', 'mut_inverted_masks', 'tests' # inv mask in top level
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
                raise ValueError('Error reading stats yaml file\n\n{}'.format(e))
    except IOError as e:
        raise IOError("cannot find or open stats config file: {}".format(config_path))
    addict_config = Dict(config)

    incorrect = unkown_options(config, TOP_LEVEL_KNOWN_OPTIONS + TOP_LEVEL_KNOWN_PATHS)
    if incorrect:
        raise ValueError('incorrect option "{}" in stats yaml file.\nDo you mean {} ?'.format(*incorrect))

    try:
        config['data']
    except KeyError:
        raise ValueError("stats config file needs a 'data' entry. Are you using the correct config file?")

    else:
        # Check each data entry for correct options (not done yet
        for key, stats_job_config in config['data'].items():
            incorrect = unkown_options(stats_job_config, STATS_JOB_KNOWN_OPTIONS)
            if incorrect:
                raise ValueError('incorrect option "{}" in stats yaml file.'
                         '\nAvaible options are \n\n{}'.format(incorrect[0], '\n'.join(STATS_JOB_KNOWN_OPTIONS)))

    ip = invalid_paths(config_path, config, TOP_LEVEL_KNOWN_PATHS)
    if len(ip) > 0:
        raise IOError('File or folder Paths do not exist\n{}'.format('\n'.join(ip)))

    return addict_config


def unkown_options(config, available_options):
    for param in config:
        if param not in available_options:
            closest_match = difflib.get_close_matches(param, TOP_LEVEL_KNOWN_OPTIONS)
            return param, closest_match
    return False


def invalid_paths(config_path, config, values):
    """
    Check that paths specified in config exist relative to config
    Parameters
    ----------
    config_path
    paths

    Returns
    -------
    list
        invalid paths

    """
    invalid = []
    config_dir = dirname(config_path)
    for v in values:
        if not config.get(v):
            continue
        option = config.get(v)
        path = join(config_dir, option)
        if not os.path.exists(path):
            invalid.append(path)
    return invalid
