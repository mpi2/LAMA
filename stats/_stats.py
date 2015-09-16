import numpy as np
import SimpleITK as sitk
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import rpy2.robjects as robj
import scipy.stats.mstats as mstats # Stats module that works on masked numpy arrays

rstats = importr('stats')


class AbstractStatisticalTest(object):
    """
    Generates the statistics. Can be all against all or each mutant against all wildtypes
    """
    def __init__(self, wt_data, mut_data, mask=None):
        """
        Parameters
        ----------
        wt_data: list
            list of ndarrays
        mask: numpy nd array
            mask array 3D
        """
        self.wt_data = wt_data
        self.mut_data = mut_data
        self.mask = mask
        self.masked_wt_data = self._get_masked_data(self.wt_data, mask)
        self.masked_mut_data = self._get_masked_data(self.mut_data, mask)

        self.filtered_tscores = False  # The final result will be stored here

    def _get_masked_data(self, data, mask=None):
        """
        Mask the numpy arrays. Numpy masked arrays can be used in scipy stats tests
        http://docs.scipy.org/doc/scipy/reference/stats.mstats.html

        If no mask, we do not mask. For eaxmple GLCM data is premasked during generation

        Parameters
        ----------
        data: list
            list of numpy 3D arrays
        Returns
        -------
        masked ndarray of arrays
        """

        flat_data = self._flatten(data)
        if mask != None:
            mask_ = 1 - mask.flatten()  # Copy the mask to get the same number of data arrays. Better way to fo this?
            flat_data = [np.ma.masked_array(a, mask_) for a in flat_data]
        return np.array(flat_data)

    @staticmethod
    def _flatten(arrays):
        one_d = []
        for arr in arrays:
            f = arr.flatten()
            one_d.append(f)
        stacked = np.vstack(one_d)
        return stacked

    def run(self):
        raise NotImplementedError


    @staticmethod
    def _result_cutoff_filter(tstats, qvalues):
        """
        Convert to numpy arrays and set to zero any tscore that has a corresponding pvalue > 0.05

        Parameters
        ----------

        """
        assert len(tstats) == len(qvalues)
        t = np.array(tstats)
        q = np.array(qvalues)
        mask = q > 0.05  # remove hard coding
        t[mask] = 0

        return t

    def get_result_array(self):
        return self.filtered_tscores



class TTest(AbstractStatisticalTest):
    """
    Compare all the mutants against all the wild type. Generate a stats overlay
    """
    def __init__(self, *args):
        super(TTest, self).__init__(*args)
        self.stats_method_object = None
        self.fdr_class = BenjaminiHochberg

    def run(self):
        """
        Returns
        -------
        sitk image:
            the stats overlay
        """

        # The masked pvalues will be Nan
        # The masked tsatistics come out as 1.0
        tstats, pvalues = mstats.ttest_ind(self.masked_mut_data, self.masked_wt_data)

        fdr = self.fdr_class(pvalues)
        qvalues = fdr.get_qvalues()

        self.filtered_tscores = self._result_cutoff_filter(tstats, qvalues)


class AbstractFalseDiscoveryCorrection(object):
    """
    Given a set of pvalues or other statistical measure, correct based on a method defined in the subclass
    """
    def __init__(self, masked_pvalues):
        """
        Parameters
        ----------
        pvalues: array
            list of pvalues to correct
        mask: numpy 3D array
        """
        self.pvalues = masked_pvalues

    def get_qvalues(self):
        raise NotImplementedError


class BenjaminiHochberg(AbstractFalseDiscoveryCorrection):
    def __init__(self, *args):
        super(BenjaminiHochberg, self).__init__(*args)

    def get_qvalues(self):
        qvals = np.array(rstats.p_adjust(FloatVector(self.pvalues), method='BH'))
        qvals[np.isnan(qvals)] = 1
        qvals[np.isneginf(qvals)] = 1
        qvals[np.isinf(qvals)] = 1
        return qvals
