"""
Test cases for the reg_stats module rewrite.
Some of these tests may well be uneeded, but I'm just trying to learn TDD here
"""

import unittest
import tempfile
from os.path import join

import yaml

from stats import bu_lama_stats


# Example part of a config file needed for the stats module
config = {'normalised_output': 'normalised_output',
          'output_dir': 'output',
          'deformations': 'deformations',
          'jacobians': 'jacobians',
          'glcms': 'glcms',
          'inverted_transforms': 'inverted_transforms',
          'inverted_isosurfaces': 'inverted_isosurfaces',
          'inverted_labels': 'inverted_labels',
          'stats': 'stats'}


class Stats(unittest.TestCase):

    def setUp(self):
        self.config_dir = tempfile.gettempdir()
        self.testfile = join(self.config_dir, 'test_statsconfig.yaml')

        with open(self.testfile, 'w+') as fh:
            fh.write(yaml.dump(config))
        self.stats = bu_lama_stats.Stats(self.testfile)


    def test_stats_loads_yaml(self):
        self.assertDictEqual(config, self.stats.get_config())

    def test_stats_factory(self):
        intensity_stats = bu_lama_stats.IntensityStats(config, self.config_dir)

        self.assertIsInstance(intensity_stats.data_getter, bu_lama_stats.ScalarDataGetter)




