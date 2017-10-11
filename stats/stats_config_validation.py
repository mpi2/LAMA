import yaml
import logging
import sys
from os.path import join, dirname
import difflib

sys.path.insert(0, join(dirname(__file__), '..'))
from lib.addict import Dict

KNOWN_OPTIONS = (
	'data', 'fixed_mask', 'label_map_path', 'littermate_controls', 'mut_staging_file', 'wt_staging_file', 'n1',
	'voxel_size', 'invert_config_file', 'project_name', 'organ_names', 'blur_fwhm', 'formulas', 'Use_auto_staging',
	'littermate_pattern'
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

	incorrect = unkown_options(config)
	if incorrect:
		sys.exit('incorrect option "{}" in stats yaml file.\nDo you mean {} ?'.format(*incorrect))

	try:
		config['data']
	except KeyError:
		logging.error("stats config file needs a 'data' entry. Are you usin gthe correct config file?")

	else:
		# Check each data entry for correct options (not done yet
		pass
	return addict_config


def unkown_options(config):
	for param in config:
		if param not in KNOWN_OPTIONS:
			closest_match = difflib.get_close_matches(param, KNOWN_OPTIONS)
			return param, closest_match
	return False
