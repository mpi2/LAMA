"""
Test cases for the reg_stats module rewrite.
Some of these tests may well be uneeded, but I'm just trying to learn TDD here
"""

import unittest
import tempfile
import yaml
from os.path import join, relpath

import reg_stats_new

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
        self.stats = reg_stats_new.Stats(self.testfile)


    def test_stats_loads_yaml(self):
        self.assertDictEqual(config, self.stats.get_config())



class StatsFactory(unittest.TestCase):
    def test_stats_factory(self):
        glcm_stats = reg_stats_new.StatsFactory('glcm')
        deformations_stats = reg_stats_new.StatsFactory('intensity')
        jacobians_stats = reg_stats_new.StatsFactory('jacobians')

        self.assertIsInstance()

