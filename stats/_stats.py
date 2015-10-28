import numpy as np
import scipy.stats.stats as stats
import gc
from os.path import join
import os.path
import tempfile
import subprocess


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

        # The masked pvalues will be Nan
        # The masked tsatistics come out as 1.0
        tstats, pvalues = self.runttest()

        # np.save('test_tstat', tstats)
        # np.save('pvale_test', pvalues.data)

        #Get the mask for the frd calculation
        # mask = np.ma.getmask(self.masked_mut_data[0])


        fdr = self.fdr_class(pvalues)
        qvalues = fdr.get_qvalues(self.mask)
        gc.collect()

        # np.save('qvalues_test', qvalues)

        #self.filtered_tscores = self._result_cutoff_filter(tstats, qvalues)
        self.filtered_tscores = self._result_cutoff_filter(tstats, qvalues) # modifies tsats in-place

    #@profile
    def runttest(self):  # seperate method for profiling
         return stats.ttest_ind(self.mut_data, self.wt_data)


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
        self.pvalues[mask == False] = 'nan'
        pval_file_for_R = join(tempfile.gettempdir(), 'pvals_for_R_temp.csv')
        qval_results_file = join(tempfile.gettempdir(), 'qvals_from_R.csv')

        np.savetxt(pval_file_for_R, self.pvalues)
        gc.collect()
        stats_dir = os.path.dirname(os.path.realpath(__file__))
        p_adjust_script = join(stats_dir, PADJUST_SCRIPT )

        # Run the Rscript
        try:
            subprocess.check_call(['Rscript', p_adjust_script, pval_file_for_R, qval_results_file])
        except subprocess.CalledProcessError:
            print "The FDR calculation failed. Is R installed and do you have 'Rscript' in your path?"
            raise

        # read in the qvals

        qvals = np.genfromtxt(qval_results_file, dtype=np.float16)  # Add dtype info


        #qvals = np.array(rstats.p_adjust(FloatVector(self.pvalues), method='BH'))
        #gc.collect()
        #self.pvalues[mask == False] = np.NaN
        #qvals = multitest(self.pvalues, method='fdr_bh')[1]
        qvals[np.isnan(qvals)] = 1
        qvals[np.isneginf(qvals)] = 1
        qvals[np.isinf(qvals)] = 1
        return qvals


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
