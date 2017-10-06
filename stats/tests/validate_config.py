"""
Reads the lama_stats yaml config into a dict.
Validates the config and return the dict if all OK
"""

import logging

def load_config_from_file(config_path):
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
    try:
        config['data']
    except KeyError:
        logging.error("stats config file needs a 'data' entry. Are you usin gthe correct config file?")

    return addict_config