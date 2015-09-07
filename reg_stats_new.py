import yaml
from os.path import relpath, join, dirname

import SimpleITK as sitk
import numpy as np

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

    def path_generator(self, path):
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
            self.stats_objects.append(IntensityStats(self.config, self.config_dir))


class AbstractPhenotypeStatistics(object):
    """
    The base class for the statistics generators
    """
    def __init__(self, config, config_dir):
        self.config = config
        self.config_dir = config_dir
        self.data_getter = self.set_data_getter()

    def set_data_getter(self):
        raise NotImplementedError

    def get_data_dir(self):
        wt_data_dir = join(self.config_dir, self.wt_data_basename)
        mut_data_dir = join(self.config, self.mut_data_basename)

class IntensityStats(AbstractPhenotypeStatistics):
    def __init__(self, *args):
        super(IntensityStats, self).__init__(*args)

    def set_data_getter(self):
        wt_data_dir = []

        return ScalarDataGetter()


class GlcmStats(AbstractPhenotypeStatistics):
    def __init__(self):
        pass


class JacobianStats(AbstractPhenotypeStatistics):
    def __init__(self):
        pass


class DeformationStats(AbstractPhenotypeStatistics):
    def __init__(self):
        pass


class AbstractDataGetter(object):
    """
    Gets the data. Could be scalar, vector or texture
    """
    def __init__(self, wt_data_dir, mut_data_dir):
        """
        Parameters
        ----------
        wt_mut_data_dir: str
            path to folder containing data volumes
        mut_data_dir: str
            path to folder containing data volumes

        Notes
        _____
        Data can exist in sub-folders
        """
        self.wt_data_dir = wt_data_dir
        self.mut_data_dir = mut_data_dir
        self.wt_paths, self.mut_paths = self.get_data_paths()
        self.wt_data, self.mut_data = self.generate_data()

    def get_data_paths(self):
        """
        Get paths to the data
        """
        wt_paths = common.GetFilePaths(self.wt_data_dir)
        mut_paths = common.GetFilePaths(self.mut_data_dir)
        return wt_paths, mut_paths

    def generate_data(self):
        """
        Process the data so it's in scalar form to be analysed by the stats classes
        """
        wt_data = []
        mut_data = []
        for wt_data_path in self.wt_paths():
            wt_data.append(self.get_data(wt_data_path))
        for mut_data_path in self.mut_paths:
            mut_data.append(self.get_data(mut_data_path))

        return wt_data, mut_data

    def get_data(self, path_):
        raise NotImplementedError


class ScalarDataGetter(AbstractDataGetter):
    """
    Processess the image intensity data:
        - The normnalised registered images
        - Spatial jacobians
    """
    def __init__(self, *args):
        super(ScalarDataGetter, self).__init__(*args)

    def get_data(self, path_):
        """
        For Intensity and jacobians, we already have scalr data that can go into the stats calculations as-is
        So just return it
        """
        return sitk.GetArrayFromImage(sitk.ReadImage(path_))

class DeformationDataGetter(AbstractDataGetter):
    """
    Process the deformations fields generated during registration
    """
    def __init__(self, *args):
        super(DeformationDataGetter, self).__init__(*args)

    def get_data(self, path_):
        """
        Calculates the deformation vector magnitude at each voxel position
        """
        arr = sitk.GetArrayFromImage(sitk.ReadImage(path_))
        return np.sqrt((arr*arr).sum(axis=3))

class GlcmDataGetter(AbstractDataGetter):
    """
    Get data from the grey level co-occurence matrices
    """
    def __init__(self, *args):
        super(GlcmDataGetter, self).__init__(*args)

    def get_data(self, path_):
        """
        The data will be in this form:
        A numpy npz for both wt and mutant. Aeach containing glcms for each block for each specimen
        """
        pass # Need to call the glcm module and extract the features

class AbstractStatsGenerator(object):
    """
    Generates the statistics. Can be all against all or each mutant against all wildtypes
    """
    def __init__(self):
        pass

class ManyAgainstManyStats(AbstractStatsGenerator):
    """
    Compare all the mutants against all the wild type. Generate a stats overlay
    """
    def __init__(self):
        pass

class OneAgainstMany(AbstractStatsGenerator):
    """
    Compare each mutant against all the wildtypes. Hopefully catches some variable penetrance
    """
    def __init__(self):
        pass
