"""
Read the stats config and do some validation.


TODO: Need an option to skip inversions
TODO: Add default options for some parameters
"""
import numbers
from pathlib import Path
from addict import Dict


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
    def path(path):
        return True if Path(path).is_file() else False

    def seq(it, allowed):
        if hasattr(it, '__iter__'):
            given = set(it)
            if len(given) != len(given.intersection(allowed)):
                return False
            else:
                return True
        else:
            return False

    def num(n, min=None, max=None):
        if not isinstance(n, numbers.Number):
            raise ValueError(f'{key} must be a number')
        wrong = False
        if min is not None:
            if n < min:
                wrong = True
        if max is not None:
            if n > max:
                wrong = True
        if wrong:
            raise ValueError(f'{key} should be a number with min {min} and max {max}')

    def bool_(b):
        return isinstance(b, bool)

    def options(given, allowed):
        if given not in allowed:
            raise ValueError(f'{given} is not an allowed option')

    schema = {
        'stats_types': {
            'required': True,
            'validate': (seq, ['intensity', 'jacobians', 'organ_volume'])
        },
        'blur_fwhm': {
            'required': False,
            'validate': (num, 0)
        },
        'voxel_size':{
            'required': False,
            'validate': (num, 0)
        },
        'mask': {
            'required': True,
            'validate': [path]
        },
        'label_info':{
            'required': False,
            'validate': [path]
        },
        'label_map':{
            'required': True,
            'validate': [path]
        },
        'reg_folder': {
            'required': False,
            'validate': [lambda x: isinstance(x, str)]
        },
        'jac_folder': {
            'required': False,
            'validate': [lambda x: isinstance(x, str)]
        },
        'invert_stats': {
            'required': False,
            'validate': [bool_]
        },
        'normalise': {
            'required': False
        },
        'use_staging': {
            'required': False,
            'validate': [options, [True, False]]  # Maybe add option to use organ volume or crown to rump length
        },
        'baseline_ids': {
            'required': False,
            'validate': [path]
        },
        'mutant_ids': {
            'required': False,
            'validate': [path]
        },
        'normalise_organ_vol_to_mask': {
            'required': False,
            'validate' : [bool_]
        }


    }

    # Check for required keys in config
    for key, data in schema.items():
        if data['required']:
            if key not in config.keys():
                raise ValueError(f'Required key "{key}" not present in config')

    # Validate the data in the config
    for key, data in config.items():

        if key not in schema:
            raise ValueError(f'{key} is not a valid stats entry\nValid keys are {schema.keys()}')

        # Do validation on values
        v = schema[key].get('validate')
        if v:
            if len(v) == 1:
                v[0](data)
            else:
                v[0](data, v[1])
