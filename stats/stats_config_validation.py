import yaml
import logging
import sys
from os.path import join, dirname
import os
import difflib
from enum import Enum

sys.path.insert(0, join(dirname(__file__), '..'))
from lib.addict import Dict


class TopLevelOptions(Enum):
    """
    Enum containing allowable options in the top level (not nested) part of the stast config
    This allows for validation and creation in phenodetect to work together as well as allowing for the chage of
    the entry value in the config to change without chaging the code
    """
    data = 'data'
    littermates = 'littermate_controls'
    do_n1 = 'n1'
    voxel_size = 'voxel_size'
    project_name = 'project_name'
    blur_fwhm = 'blur_fwhm'
    formulas = 'formulas'
    use_auto_staging = 'use_auto_staging'
    littermate_pattern = 'littermate_pattern'
    mutant_ids = 'mutant_ids'
    wildtype_ids = 'wildtype_ids'


class TopLevelPaths(Enum):
    output_dir='output_dir'
    fixed_mask = 'fixed_mask'
    label_map = 'label_map'
    mut_staging_file = 'mut_staging_file'
    wt_staging_file = 'wt_staging_file'
    invert_config_file = 'invert_config_file'
    label_names = 'label_names'
    inverted_masks = 'inverted_masks'


class StatsEntryOptions(Enum):
    wt = 'wt'
    mut = 'mut'
    normalisation = 'normalisation'
    wt_inverted_masks = 'wt_inverted_masks'
    mut_inverted_masks = 'mut_inverted_masks'
    tests = 'tests'



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

    incorrect = unkown_options(config, [x.value for x in TopLevelOptions] + [x.value for x in TopLevelPaths])
    if incorrect:
        raise ValueError('incorrect option "{}" in stats yaml file.\nDo you mean {} ?'.format(*incorrect))

    try:
        config['data']
    except KeyError:
        raise ValueError("stats config file needs a 'data' entry. Are you using the correct config file?")

    else:
        # Check each data entry for correct options (not done yet
        for key, stats_job_config in config['data'].items():
            incorrect = unkown_options(stats_job_config, [x.value for x in StatsEntryOptions])
            if incorrect:
                raise ValueError('incorrect option "{}" in stats yaml file.'
                         '\nAvaible options are \n\n{}'.format(incorrect[0], '\n'.join([x.value for x in StatsEntryOptions])))

    ip = invalid_paths(config_path, config, TopLevelPaths)
    if len(ip) > 0:
        raise IOError('File or folder Paths do not exist\n{}'.format('\n'.join(ip)))

    return addict_config


def unkown_options(config, available_options):
    for param in config:
        if param not in available_options:
            closest_match = difflib.get_close_matches(param, available_options)
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
