import numpy as np
import scipy.stats.stats as stats
import gc
from os.path import join
import os.path
import csv
from collections import defaultdict
import tempfile
import subprocess
import sys
import struct
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import rpy2.robjects as robj
rstats = importr('stats')
# import statsmodels.stats.

MINMAX_TSCORE = 50 # If we get very large tstats or in/-inf this is our new max/min
# PADJUST_SCRIPT = 'r_padjust.R'
LINEAR_MODEL_SCIPT = 'lmFast.R'
FDR_SCRPT = 'r_padjust.R'
VOLUME_METADATA_NAME = 'volume_metadata.csv'
DATA_FILE_FOR_R_LM = 'tmp_data_for_lm'
PVAL_R_OUTFILE = 'pvals_out.dat'
GROUPS_FILE_FOR_LM = 'groups.csv'

class AbstractStatisticalTest(object):
    """
    Generates the statistics. Can be all against all or each mutant against all wildtypes
    """
    def __init__(self, wt_data, mut_data, mask, zscores, groups=None):
        """
        Parameters
        ----------
        wt_data: list
            list of masked 1D ndarrays
        mut_data: list
            list of masked 1D ndarrays
        zscores: np.ndarry (1d)
            Number of standard deviations between wildtypes and mutants at each pixel
        groups: dict/None
            For linear models et. al. contains groups membership for each volume
        """
        self.groups = groups
        self.mask = mask
        self.wt_data = wt_data
        self.mut_data = mut_data
        self.filtered_tscores = False  # The final result will be stored here
        self.zscores = zscores

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

    def get_volume_metadata(self):
        """
        Get the metada for the volumes, such as sex and (in the future) scan date
        Not currently used
        """
        def get_from_csv(csv_path):
            with open(csv_path, 'rb') as fh:
                reader = csv.reader(fh, delimiter=',')
                first = True
                for row in reader:
                    if first:
                        first = False
                        header = row
                    else:
                        vol_id = row[0]
                        for i in range(1, len(row)):
                            meta_data[vol_id][header[i]] = row[i]

        mut_vol_metadata_path = join(self.mut_proj_dir, self.in_dir, VOLUME_METADATA_NAME)
        wt_vol_metadata_path = join(self.wt_config_dir, self.wt_config['inputvolumes_dir'], VOLUME_METADATA_NAME)

        if not os.path.exists(mut_vol_metadata_path) or not os.path.exists(wt_vol_metadata_path):
            print 'Cannot find volume metadata, will only be able to do linear model anlysis with genotype'
            return False

        meta_data = defaultdict(dict)

        get_from_csv(mut_vol_metadata_path)
        get_from_csv(wt_vol_metadata_path)

        return meta_data


class LinearModelR(AbstractStatisticalTest):
    def __init__(self, *args):
        super(LinearModelR, self).__init__(*args)
        self.stats_method_object = None  # ?
        self.fdr_class = BenjaminiHochberg

    def set_groups_file(self, file_path):
        """
        Set the path to the groups csv that specifies which groups each volume belongs to.
        If this is not set, we will fall back to just using genotype as this can be determined without it.
        """

    def run(self):
        # These contain the chunked stats results
        pvals = []

        size = self.wt_data[0].size

        # np.array_split provides a split view on the array so does not increase memory
        # The result will be a bunch of arrays split across the second dimension
        groups = ['wildtype'] * len(self.wt_data)
        groups.extend(['mutant'] * len(self.mut_data))

        # Write out the groups/levels csv
        groups_file = join(tempfile.gettempdir(), GROUPS_FILE_FOR_LM)
        with open(groups_file, 'w') as gf:
            for g in groups:
                gf.write(g + '\n')

        linear_model_script = join(os.path.dirname(os.path.realpath(__file__)), LINEAR_MODEL_SCIPT)

        pval_out_file = join(tempfile.gettempdir(), PVAL_R_OUTFILE)

        data = np.vstack((self.wt_data, self.mut_data))
        data_filtered = np.hstack(data[:, np.argwhere(self.mask != False)]).T

        num_pixels = data_filtered.shape[1]
        chunk_size = 200000
        num_chunks = num_pixels / chunk_size
        print 'number of chuncks: ', num_chunks

        # Loop over the data in chunks
        chunked_data = np.array_split(data_filtered, num_chunks, axis=1)

        i = 0
        for data_chucnk in chunked_data:
            i += 1
            print 'done LM chunks', i
            pixel_file = join(tempfile.gettempdir(), DATA_FILE_FOR_R_LM)
            numpy_to_dat(np.vstack(data_chucnk), pixel_file)

            # fit the data to a linear m odel and extrat the tvalue
            try:
                subprocess.check_output(['Rscript',
                                         linear_model_script,
                                         pixel_file,
                                         groups_file,
                                         pval_out_file])
            except subprocess.CalledProcessError:
                print "r linear model failed"
                raise

            # Read in the pvalue results
            p = np.fromfile(pval_out_file, dtype=np.float64)

            # Convert all NANs in the pvalues to 1.0. Need to check that this is appropriate
            p[np.isnan(p)] = 1.0
            pvals.append(p)

        #tstats_array = np.hstack(tstats)
        pvals_array = np.hstack(pvals)

        result_pvals = np.ones(size)
        result_pvals[self.mask != False] = pvals_array

        fdr = self.fdr_class(result_pvals)
        qvalues = fdr.get_qvalues(self.mask)

        gc.collect()

        filtered_zscores = self._result_cutoff_filter(self.zscores, qvalues)

        print
        print 'filt min', filtered_zscores.min()
        print 'filt max', filtered_zscores.max()
        print 'filt mean', filtered_zscores.mean()

        # Remove infinite values
        filtered_zscores[np.isnan(filtered_zscores)] = 0.0
        filtered_zscores[filtered_zscores > MINMAX_TSCORE] = MINMAX_TSCORE
        filtered_zscores[filtered_zscores < -MINMAX_TSCORE] = - MINMAX_TSCORE

        self.filtered_tscores = filtered_zscores

        # Put the tscores back into the correct shapes overlay



class TTest(AbstractStatisticalTest):
    """
    Compare all the mutants against all the wild type. Generate a stats overlay

    When working with
    """
    def __init__(self, *args):
        super(TTest, self).__init__(*args)
        self.stats_method_object = None #?
        self.fdr_class = BenjaminiHochberg

    def run(self):
        """
        Returns
        -------
        sitk image:
            the stats overlay
        """

        # Temp. just split out csv of wildtype and mutant for R


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
            pval_chunk = pval_chunk.astype(np.float32)
            tstats.extend(tstats_chunk)
            pvals.extend(pval_chunk)

        pvals = np.array(pvals)
        tstats = np.array(tstats)

        fdr = self.fdr_class(pvals)
        qvalues = fdr.get_qvalues(self.mask)
        gc.collect()

        #self.filtered_tscores = self._result_cutoff_filter(tstats, qvalues)
        self.filtered_tscores = self._result_cutoff_filter(tstats, qvalues) # modifies tsats in-place

        # Remove infinite values
        self.filtered_tscores[self.filtered_tscores > MINMAX_TSCORE] = MINMAX_TSCORE
        self.filtered_tscores[self.filtered_tscores < -MINMAX_TSCORE] = - MINMAX_TSCORE

    #@profile
    def runttest(self, wt, mut):  # seperate method for profiling

        return stats.ttest_ind(mut, wt)

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

        # pvals_corrected[pvals_corrected > 1] = 1
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


class BenjaminiHochbergR(AbstractFalseDiscoveryCorrection):
    def __init__(self, *args):
        super(BenjaminiHochbergR, self).__init__(*args)

    def get_qvalues(self, mask):
        """
        Mask ndarray of booleans. True == masked
        """
        print 'Doing r calculation'
        self.pvalues[mask == False] = robj.NA_Real
        qvals = np.array(rstats.p_adjust(FloatVector(self.pvalues), method='BH'))
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


def numpy_to_dat(mat, outfile):

    # create a binary file
    binfile = file(outfile, 'wb')
    # and write out two integers with the row and column dimension

    header = struct.pack('2I', mat.shape[0], mat.shape[1])
    binfile.write(header)
    # then loop over columns and write each
    for i in range(mat.shape[1]):
        data = struct.pack('%id' % mat.shape[0], *mat[:, i])
        binfile.write(data)

    binfile.close()
