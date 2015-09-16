
from os.path import join
import sys
import os
# Hack. Relative package imports won't work if this module is run as __main__
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common
import SimpleITK as sitk
from invert import BatchInvertLabelMap
from _data_getters import GlcmDataGetter, DeformationDataGetter, ScalarDataGetter
import numpy as np


class AbstractPhenotypeStatistics(object):
    """
    The base class for the statistics generators
    """
    def __init__(self, config_dir, out_dir, wt_data_dir, mut_data_dir, mask_array=None):
        """
        Parameters
        ----------
        mask_array: numpy ndarray
            3D mask array
        """
        self.config_dir = config_dir
        self.out_dir = out_dir
        common.mkdir_if_not_exists(self.out_dir)
        self.mask = mask_array
        self._wt_data_dir = wt_data_dir
        self._mut_data_dir = mut_data_dir

        self.wt_data = None
        self.mut_data = None

        # These need to be set from outside
        self.tests = None
        self.shape = None

    def _set_data(self):
        self.dg = dg = self.data_getter(self._wt_data_dir, self._mut_data_dir)
        self.shape = dg.shape
        self.wt_data = dg.wt_data
        self.mut_data = dg.mut_data

    def run(self, stats_object, analysis_prefix):
        self._set_data()
        so = stats_object(self.wt_data, self.mut_data, self.mask)
        so.run()
        stats_array = so.get_result_array()
        reshaped_array = self._reshape_data(stats_array)
        result_img = sitk.GetImageFromArray(reshaped_array)
        outfile = join(self.out_dir, analysis_prefix + '.nrrd')  # remove hard coding of nrrd
        sitk.WriteImage(result_img, outfile)

    def _reshape_data(self, result_data):
        """
        The data normally can be reshaped straight from the stats analysis. In that case leave as is
        If data has been chunked or subsampled as in GLCM analysis, this method needs overriding

        Parameters
        ----------
        result_data: ndarray
            possibly a 1d array that needs reshaping
        """
        return result_data.reshape(self.original_shape)


class IntensityStats(AbstractPhenotypeStatistics):
    def __init__(self, *args):
        super(IntensityStats, self).__init__(*args)
        self.data_getter = ScalarDataGetter


class GlcmStats(AbstractPhenotypeStatistics):
    def __init__(self, *args):
        super(GlcmStats, self).__init__(*args)
        self.data_getter = GlcmDataGetter

    def _reshape_data(self, result_data):
        """
        The data from the GLCM analysis is subsampled and so smaller than the original data. To be able to overlay
        onto real image data, we need to upsample the result

        Parameters
        ----------

        Returns
        -------
        A numpy ndarray? Should be 1D
        """
        shape = self.shape
        chunk_size = self.dg.glcm_chunk_size
        out_array = np.zeros(self.shape)
        i = 0
        for z in range(0, shape[0] - chunk_size, chunk_size):
            for y in range(0, shape[1] - chunk_size, chunk_size):
                for x in range(0, shape[2] - chunk_size, chunk_size):
                    out_array[z: z + chunk_size, y: y + chunk_size, x: x + chunk_size] = result_data[i]
                    i += 1

        return out_array


class JacobianStats(AbstractPhenotypeStatistics):
    def __init__(self):
        self.data_getter = ScalarDataGetter


class DeformationStats(AbstractPhenotypeStatistics):
    def __init__(self):
        pass





