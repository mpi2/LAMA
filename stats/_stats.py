import numpy as np
import scipy.stats.stats as stats
from scikits.statsmodels.stats import multipletests as multtest
import gc
from os.path import join
import os.path
import tempfile
import subprocess

MINMAX_TSCORE = 50 # If we get very large tstats or in/-inf this is our new max/min
PADJUST_SCRIPT = 'r_padjust.R'

class AbstractStatisticalTest(object):
    """
    Generates the statistics. Can be all against all or each mutant against all wildtypes
    """
    def __init__(self, wt_data, mut_data, mask):
        """
        Parameters
        ----------
        wt_data: list
            list of masked 1D ndarrays
        mut_data: list
            list of masked 1D ndarrays
        """
        self.mask = mask
        self.wt_data = wt_data
        self.mut_data = mut_data
        self.filtered_tscores = False  # The final result will be stored here

    def run(self):
        raise NotImplementedError


    @staticmethod
    def _result_cutoff_filter(t, q):
        """
        Convert to numpy arrays and set to zero any tscore that has a corresponding pvalue > 0.05

        Parameters
        ----------

        """
        assert len(t) == len(q)
        # t = np.array(tstats)
        # q = np.array(qvalues)
        mask = q > 0.05  # remove hard coding
        t[mask] = 0

        return t

    def get_result_array(self):
        return self.filtered_tscores

class TTest(AbstractStatisticalTest):
    """
    Compare all the mutants against all the wild type. Generate a stats overlay

    When working with
    """
    def __init__(self, *args):
        super(TTest, self).__init__(*args)
        self.stats_method_object = None
        self.fdr_class = BenjaminiHochberg

    #@profile
    def run(self):
        """
        Returns
        -------
        sitk image:
            the stats overlay
        """

        # Chunk the ttest to save memory

        # These contain the chunked stats results
        tstats = []
        pvals = []



        # np.array_split provides a split view on the array so does not increase memory
        # The result will be a bunch of arrays split down the second dimension
        chunked_mut = np.array_split(self.mut_data, 10, axis=1)
        chunked_wt = np.array_split(self.wt_data, 10, axis=1)

        for wt_chunks, mut_chunks in zip(chunked_wt, chunked_mut):



            tstats_chunk, pval_chunk = self.runttest(wt_chunks, mut_chunks)
            pval_chunk[np.isnan(pval_chunk)] = 0.1
            print pval_chunk.min(), pval_chunk.mean()
            pval_chunk = pval_chunk.astype(np.float32)
            tstats.extend(tstats_chunk)
            pvals.extend(pval_chunk)

        pvals = np.array(pvals)
        tstats = np.array(tstats)


        fdr = self.fdr_class(pvals)
        qvalues = fdr.get_qvalues(self.mask)
        gc.collect()

        # np.save('qvalues_test', qvalues)

        #self.filtered_tscores = self._result_cutoff_filter(tstats, qvalues)
        self.filtered_tscores = self._result_cutoff_filter(tstats, qvalues) # modifies tsats in-place

        # Remove infinite values
        self.filtered_tscores[self.filtered_tscores > MINMAX_TSCORE] = MINMAX_TSCORE
        self.filtered_tscores[self.filtered_tscores < -MINMAX_TSCORE] = - MINMAX_TSCORE

    #@profile
    def runttest(self, wt, mut):  # seperate method for profiling

        for i in range(wt[0].size):
            w_test = []
            for w in wt:
                w_test.append(w[i])
            m_test = []
            for m in mut:
                m_test.append(m[i])
            tscore, pvalue = stats.ttest_ind(m_test, w_test)
            if tscore > 10:
                print tscore, pvalue, m_test, w_test



        #return stats.ttest_ind(mut, wt)

    def split_array(self, array):
        """
        Split array into equal-sized chunks + remainder
        """
        return np.array_split(array, 5)


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

    #@profile
    def get_qvalues(self, mask):
        """
        Mask ndarray of booleans. True == masked
        """

        # Write out pvalues to temporary file for use in R
        pvals = self.pvalues[mask != False]
        size = self.pvalues.size

        pvals_sortind = np.argsort(pvals)
        pvals_sorted = pvals[pvals_sortind]
        sortrevind = pvals_sortind.argsort()

        ecdffactor = self.ecdf(pvals_sorted)

        pvals_corrected_raw = pvals_sorted / ecdffactor


        pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]

        pvals_corrected[pvals_corrected > 1] = 1
        pvals_corrected[np.isnan(pvals_corrected)] = 1
        pvals_corrected[np.isneginf(pvals_corrected)] = 1
        pvals_corrected[np.isinf(pvals_corrected)] = 1

        pvals_resorted = pvals_corrected[sortrevind]

        # Insert the mask positions back into the array
        result = np.ones(size) ## add 16f dtype?
        result[mask != False] = pvals_resorted

        return result


    def ecdf(self, x):
        '''no frills empirical cdf used in fdrcorrection
        '''
        nobs = len(x)
        return np.arange(1,nobs+1)/float(nobs)


class OneAgainstManytest(object):
    def __init__(self, wt_data, zscore_cutoff=3):
        """
        Perform a pixel-wise z-score analysis of mutants compared to a set of  wild types

        Parameters
        ----------
        wt_data: list(np.ndarray)
            list of 1d wt data
        """
        self.wt_data = wt_data
        self.zscore_cutoff = zscore_cutoff

    def process_mutant(self, mut_data):
        """
        Get the pixel-wise z-score of a mutant

        Parameters
        ----------
        mut_data: numpy ndarray
            1D masked array

        Returns
        -------
        1D np.ndarray of zscore values

        """

        z_scores = stats.zmap(mut_data, self.wt_data)

        # Filter out any values below x standard Deviations
        z_scores[np.absolute(z_scores) < self.zscore_cutoff] = 0

        # Remove nans
        z_scores[np.isnan(z_scores)] = 0

        return z_scores
