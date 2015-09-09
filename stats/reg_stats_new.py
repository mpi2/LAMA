import yaml
from os.path import relpath, join, dirname

import SimpleITK as sitk
import numpy as np
import sys
import os

from phenotype_statistics import DeformationStats, GlcmStats, IntensityStats

# Hack. Relative package imports won't work if this module is run as __main__
sys.path.insert(0, os.path.abspath('..'))
import common


class Stats(object):
    """
    Takes a stats.yaml config file and creates apropriate PhenotypeStatitics subclasses based on which analysis is to
    be performed
    """
    def __init__(self, config_path):
        self.config_path = config_path
        self.config = self.get_config()
        self.config_dir = dirname(self.config_path)
        self.stats_objects = []
        self.outdir = self.make_path(self.config['stats'])
        common.mkdir_if_not_exists(self.outdir)
        self.mask_path = self.make_path(self.config['fixed_mask'])

    def make_path(self, path):
        """
        All paths are relative to the config file dir.
        Return relative paths to the config dir
        """
        return relpath(path, self.config_dir)

    def get_config(self):
        with open(self.config_path) as fh:
            config = yaml.load(fh)
        return config

    def make_stats_objects(self):
        """
        Build the regquired stats classes for each data type
        """

        if self.config['normalised_output']:
            IntensityStats(self.config['registered_normalised'], self.config_dir, self.outdir, self.mask_path)


