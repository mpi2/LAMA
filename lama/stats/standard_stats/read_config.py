"""
Read the stats config and do some validation.
"""

from pathlib import Path
from addict import Dict
import yaml
# from strictyaml import load, Map, Str, Int, Seq, Optional, YAMLError


def read(config_path: Path) -> Dict:

    with open(config_path, 'r') as fh:

        config_yaml = yaml.load(fh)
        config = Dict(config_yaml)

        validate(config)
        return config


def validate(config: Dict):
    """
    Check for valid entries in the stats config

    Parameters
    ----------
    config
        The stats config

    Raises
    ------
    ValueError if the config file in invalid

    """

    #'allowed': ['intensity', 'jacobians', 'organ_volume']

    schema = {'stats_types': {
        'type_': list,
        'required': True,
        'allowed': ['intensity', 'jacobians', 'organ_volume']
    }
    }

    for key, data in config.items():

        if key not in schema:
            raise ValueError(f'{key} is anot a valid stats entry')

        if not isinstance(data, schema[key]['type_']):
            raise ValueError(f'The key {key} needs to map to type {schema[key]["type_"]}')

        # Check that all values in the list are in the allowed values
        if schema[key]['type_'] == list:
            given = set(data)
            allowed = set(schema[key]['allowed'])
            if len(given) != len(given.intersection(allowed)):
                raise ValueError(f'{key} maps to a list that can contain only {allowed}')
