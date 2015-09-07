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
            self.stats_objects.append(IntensityStats(self.config['registered_normalised'],
                                                     self.config_dir, self.outdir, self.mask_path))


class AbstractPhenotypeStatistics(object):
    """
    The base class for the statistics generators
    """
    def __init__(self, config, config_dir, outdir, mask_path):
        self.config = config
        self.config_dir = config_dir
        self.data_getter = self.set_data_getter()
        self.outdir = outdir
        self.mask = common.img_path_to_array(mask_path)

    def run(self):
        many_against_many = ManyAgainstManyStats(self.data_getter.wt_data, self.data_getter.mut_data, self.mask)
        many_against_many.run()

        # Todo: one against one analysis

    def set_data_getter(self):
        raise NotImplementedError


class IntensityStats(AbstractPhenotypeStatistics):
    def __init__(self, *args):
        super(IntensityStats, self).__init__(*args)

    def set_data_getter(self):
        wt_data_dir = join(self.config_dir, self.config['wt'])
        mut_data_dir = join(self.config, self.config['mut'])
        return ScalarDataGetter(wt_data_dir, mut_data_dir)


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

    @staticmethod
    def blur_volume(img):
        """
        Create memory-mapped volumes
        """
        # previous: 1.0, 8, 0.001
        blurred = sitk.DiscreteGaussian(img, 0.5, 4, 0.01, False)
        return blurred


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
        return self.blur(sitk.GetArrayFromImage(sitk.ReadImage(path_)))


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
    def __init__(self, wt_data, mut_data, mask):
        self.wt_data = wt_data
        self.mut_data = mut_data
        self.mask = mask
        self.shape = self.set_shape()

    def run(self):
        raise NotImplementedError

    def set_shape(self):
        return self.wt_data[0].shape[0:3]


class ManyAgainstManyStats(AbstractStatsGenerator):
    """
    Compare all the mutants against all the wild type. Generate a stats overlay
    """
    def __init__(self, *args):
        super(ManyAgainstManyStats, self).__init__(*args)

    def run(self):

        tstats, pvalues = ttest(self.wt_data, self.mut_data)

        # print("Calculating FDR")
        qvalues = fdr(pvalues, mask_arr)

        # print 'q', min(qvalues), max(qvalues), np.mean(qvalues)
        # reshape

        filtered_t = filter_tsats(tstats, qvalues)

        t_vol = filtered_t.reshape(shape)

        t_img = sitk.GetImageFromArray(t_vol)
        outfile = os.path.join(analysis_dir, TSCORE_OUT_SUFFIX)

        sitk.WriteImage(t_img, outfile)

    def ttest(self, wt, mut):
        """
        Calculate the pvalue and the tstatistic for the wt and mut subsampled region

        Args:
           wt (list):  wildtype values
           mut (list): mutant values

        Returns:
           tuple: (pvalue(float), is_tscore_positive(bool)
        """
        wt_flat = flatten(wt)
        mut_flat = flatten(mut)

        tscores, pvals = stats.ttest_ind(mut_flat, wt_flat)

        mask_t = np.isnan(tscores)
        tscores[mask_t] = 0
        mask_p = np.isnan(pvals)
        pvals[mask_p] = 1

        return tscores, pvals

class OneAgainstMany(AbstractStatsGenerator):
    """
    Compare each mutant against all the wildtypes. Hopefully catches some variable penetrance
    """
    def __init__(self):
        pass

