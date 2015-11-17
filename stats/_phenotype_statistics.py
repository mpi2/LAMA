
from os.path import join, basename
import sys
import os
# Hack. Relative package imports won't work if this module is run as __main__
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common
import SimpleITK as sitk
from invert import InvertVol
from _stats import OneAgainstManytest
from _data_getters import GlcmDataGetter, DeformationDataGetter, ScalarDataGetter, JacobianDataGetter
import numpy as np
import gc
import csv

STATS_FILE_SUFFIX = '_stats_'

class AbstractPhenotypeStatistics(object):
    """
    The base class for the statistics generators
    """
    def __init__(self, config_dir, out_dir, wt_data_dir, mut_data_dir, project_name, mask_array=None, groups=None, formulas=None):
        """
        Parameters
        ----------
        mask_array: numpy ndarray
            1D mask array
        groups: dict
            specifies which groups the data volumes belong to (for linear model etc.)
        """
        self.project_name = project_name
        self.config_dir = config_dir
        self.out_dir = out_dir
        common.mkdir_if_not_exists(self.out_dir)
        self.mask = mask_array  # this is a flat binary array
        self.groups = groups
        self.formulas = formulas
        self._wt_data_dir = wt_data_dir
        self._mut_data_dir = mut_data_dir

        # Obtained from the datagetter
        self.shape = None

        self.n1_stats_output = []  # Paths to the n1 anlaysis output. Use din inverting stats volumes

    def _set_data(self):
        """
        Set the wt and mut data. What are the types?
        """
        if self.groups:
            vol_order = self.get_volume_order()
        self.dg = dg = self.data_getter(self._wt_data_dir, self._mut_data_dir, self.mask, vol_order)
        self.shape = dg.shape

    def get_volume_order(self):
        """

        Returns
        -------

        list: order of volumes in groups file
        """
        if self.groups:
            order = []
            with open(self.groups, 'r') as fh:
                first = True
                reader = csv.reader(fh)
                for row in reader:
                    if first:
                        first = False
                        continue
                    else:
                        order.append(row[0])
            return order
        else:
            return None

    def run(self, stats_object, analysis_prefix):
        self._set_data()
        self._many_against_many(stats_object, analysis_prefix)
        self._one_against_many(analysis_prefix)
        del self.dg
        gc.collect()

    def _one_against_many(self, analysis_prefix):
        """
        Compare each mutant seperatley against all wildtypes
        """
        n1 = OneAgainstManytest(self.dg.wt_data)
        out_dir = join(self.out_dir, 'n1')
        common.mkdir_if_not_exists(out_dir)

        for path, mut_data in zip(self.dg.mut_paths, self.dg.mut_data):
            result = n1.process_mutant(mut_data)
            reshaped_data = self._reshape_data(result)
            out_path = join(out_dir, analysis_prefix + STATS_FILE_SUFFIX + os.path.basename(path))
            self.n1_stats_output.append(out_path)
            outimg = sitk.GetImageFromArray(reshaped_data)
            sitk.WriteImage(outimg, out_path, True)  # Compress output
        del n1
        gc.collect()

    def _many_against_many(self, stats_object, analysis_prefix):
        """
        Compare all mutants against all wild types
        """
        stats_prefix = self.project_name + '_' + analysis_prefix
        so = stats_object(self.dg.wt_data, self.dg.mut_data, self.mask,
                          self.dg.zscore_overlay, self.shape, self.out_dir, stats_prefix, self.groups, self.formulas)
        so.run()

        #stats_array = so.get_result_array()

        # # Al this to go into the stats objewct
        # reshaped_array = self._reshape_data(stats_array)
        # result_img = sitk.GetImageFromArray(reshaped_array)
        # outfile = join(self.out_dir, analysis_prefix + STATS_FILE_SUFFIX + '.nrrd')  # remove hard coding of nrrd
        # sitk.WriteImage(result_img, outfile, True)  # Stats output compresses well
        del so
        gc.collect()

    def invert(self, invert_config_path):
        """
        Invert the stats back onto the rigidly aligned volumes

        Parameters
        ----------
        invert_order: dict
            Contains inversion order information
        """
        invert_out_dir = join(self.out_dir, 'inverted')
        common.mkdir_if_not_exists(invert_out_dir)
        for stats_vol_path in self.n1_stats_output:
            single_invert_out = join(invert_out_dir, basename(stats_vol_path))
            InvertVol(invert_config_path, stats_vol_path, single_invert_out)

    def _reshape_data(self, result_data):
        """
        The data normally can be reshaped straight from the stats analysis. In that case leave as is
        If data has been chunked or subsampled as in GLCM analysis, this method needs overriding

        Parameters
        ----------
        result_data: ndarray
            possibly a 1d array that needs reshaping
        """
        return result_data.reshape(self.shape)


class IntensityStats(AbstractPhenotypeStatistics):
    def __init__(self, *args):
        super(IntensityStats, self).__init__(*args)
        self.data_getter = ScalarDataGetter


class GlcmStats(AbstractPhenotypeStatistics):
    def __init__(self, *args):
        super(GlcmStats, self).__init__(*args)
        self.data_getter = GlcmDataGetter # Currently just gets contrast measure

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

    def _mask_data(self, data):
        """
        Mask the numpy arrays. Numpy masked arrays can be used in scipy stats tests
        http://docs.scipy.org/doc/scipy/reference/stats.mstats.html

        If no mask, we do not mask. For eaxmple GLCM data is premasked during generation?????

        Parameters
        ----------
        data: list
            list of numpy 3D arrays
        Returns
        -------
        masked 1D ndarray of arrays
        """

        masked_data = np.ma.masked_where(data, np.isnan(data))
        return masked_data


class JacobianStats(AbstractPhenotypeStatistics):
    # Not used. Intensity and jacoabian analysis is the same
    def __init__(self, *args):
        super(JacobianStats, self).__init__(*args)
        self.data_getter = JacobianDataGetter


class DeformationStats(AbstractPhenotypeStatistics):
    def __init__(self, *args):
        super(DeformationStats, self).__init__(*args)
        self.data_getter = DeformationDataGetter





