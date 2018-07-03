import pandas as pd
import scipy.stats.stats as scipystats
from scipy.special import stdtr
from scipy.stats import circmean, circvar, circstd, ttest_ind, norm
import gc
from os.path import join
import os.path
import numpy as np
from scipy.stats import t as t_
import subprocess
import sys
import struct
import logging
import tempfile
from collections import defaultdict
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common
import statsmodels.stats.multitest as multitest
MINMAX_TSCORE = 50 # If we get very large tstats or in/-inf this is our new max/min
PADJUST_SCRIPT = 'rscripts/r_padjust.R'
#PADJUST_SCRIPT = 'rscripts/r_qvalues.R'
LINEAR_MODEL_SCRIPT = 'rscripts/lmFast.R'
CIRCULAR_SCRIPT = 'circular.R'
VOLUME_METADATA_NAME = 'volume_metadata.csv'
DATA_FILE_FOR_R_LM = 'tmp_data_for_lm'
PVAL_R_OUTFILE = 'tmp_pvals_out.dat'
TVAL_R_OUTFILE = 'tmp_tvals_out.dat'
GROUPS_FILE_FOR_LM = 'groups.csv'
STATS_FILE_SUFFIX = '_stats_'
PVAL_DIST_IMG_FILE = 'pval_distribution.png'
R_CHUNK_SIZE = 2000000


class AbstractStatisticalTest(object):
    """
    Generates the statistics. Can be all against all or each mutant against all wildtypes
    """
    def __init__(self, wt_data, mut_data, shape, outdir):
        """
        Parameters
        ----------
        wt_data: list
            list of masked 1D ndarrays
        mut_data: list
            list of masked 1D ndarrays
        shape: tuple
            The shape of the final result of the stats (z,y,x)
        groups: dict/None
            For linear models et. al. contains groups membership for each volume
        """
        self.outdir = outdir
        self.shape = shape
        self.wt_data = wt_data
        self.mut_data = mut_data
        self.filtered_tscores = False  # The final result will be stored here

    def run(self):
        raise NotImplementedError

    def get_result_array(self):
        return self.filtered_tscores


    def write_result(self, result_array, outpath):
        """
        """
         # Create a full size output array
        size = np.prod(self.shape)
        full_output = np.zeros(size)

        # Insert the result p and t vals back into full size array
        full_output[self.mask != False] = result_array

        reshaped_results = full_output.reshape(self.shape)
        common.write_array(reshaped_results, outpath)


class LinearModelNumpy(AbstractStatisticalTest):

    def __init__(self, *args):
        super(LinearModelNumpy, self).__init__(*args)
        self.fdr_class = BenjaminiHochberg
        self.STATS_NAME = 'LinearModelPython'

        # The output of the tests
        self.tstats = None
        self.qvals = None
        self.fdr_tstats = None

    def run(self):
        data = np.vstack((self.wt_data, self.mut_data))

        from sklearn.linear_model import LinearRegression
        n = len(data)
        genotype = ([0] * len(self.wt_data) + [1] * len(self.mut_data))
        genotype = np.array(genotype).reshape((data.shape[0], 1))

        x = genotype

        y = data

        model = LinearRegression(fit_intercept=True, copy_X=True)
        fit = model.fit(x, y)
        # So we have 2 coeffiencients for each sample first is for genotype second is for crl
        pred = fit.predict(genotype)

        # Simple linear regression
        # Test until I do maean on column
        mean_x = np.mean(genotype)

        se_slope = np.sqrt(np.sum((y - pred) ** 2, axis=0) / (n - 2)) / np.sqrt(
            np.sum(((genotype - mean_x) ** 2), axis=0))

        coef = fit.coef_.flatten()
        t = coef / se_slope
        # p = t_.cdf(t, n - 2) *2  # *2 for two sided tests
        p = t_.sf(np.abs(t), n - 2) * 2

        self.tstats = t
        self.pvals = p
        fdr_runner = self.fdr_class(p)
        self.qvals = fdr_runner.get_qvalues()

    def set_formula(self, formula):
        pass

    def set_groups(self, _):
        pass


class StatsTestR(AbstractStatisticalTest):

    def __init__(self, *args):
        super(StatsTestR, self).__init__(*args)
        self.fdr_class = BenjaminiHochberg
        self.tstats = None
        self.qvals = None

    def set_formula(self, formula):
        self.formula = formula

    def set_groups(self, groups):
        self.groups = groups

    def run(self):

        if not self.groups:
            # We need groups file for linera model
            logging.warning('linear model failed. We need groups file')
            return

        # np.array_split provides a split view on the array so does not increase memory
        # The result will be a bunch of arrays split across the second dimension

        # Make temp files for storing the output of the line levl LM results
        line_level_pval_out_file = tempfile.NamedTemporaryFile().name
        line_level_tstat_out_file = tempfile.NamedTemporaryFile().name

        # Make temp files for each mutant to store the specimen-level t-statistics
        # Not getting p-values as they ar enot needed, aI can calculate a t-cuttoff and threshold using that

        data = np.vstack((self.wt_data, self.mut_data))

        num_pixels = data.shape[1]

        num_chunks = num_pixels / R_CHUNK_SIZE
        if num_pixels < R_CHUNK_SIZE:
            num_chunks = 1
        logging.info('using formula {}'.format(self.formula))
        print 'num chunks', num_chunks

        # Loop over the data in chunks
        chunked_data = np.array_split(data, num_chunks, axis=1)
        num_data_points = chunked_data[0].shape[1]  # TThe number of voxel position or labels in a given chunk

        # These contain the chunked stats results
        line_level_pvals = []
        line_level_tvals = []

        # These contain the specimen-level statstics
        specimen_pvals = defaultdict(list)
        specimen_tstats = defaultdict(list)

        # Load in the groups file so we can get the specimne names
        groups_df = pd.read_csv(self.groups)
        mutants_df = groups_df[groups_df.genotype == 'wildtype']

        # Loop over the mutants and make temporary file for them
        specimen_level_temp_files = {}
        for i, row in mutants_df.iterrows():
            specimen_id = row['volume_id']
            specimen_level_temp_files[specimen_id] = tempfile.NamedTemporaryFile().name

        i = 0
        voxel_file = tempfile.NamedTemporaryFile().name
        for data_chucnk in chunked_data:
            print 'chunk: {}'.format(i)
            i += 1

            #stacked = np.vstack(data_chucnk)
            numpy_to_dat(data_chucnk, voxel_file)
            # import shutil
            # shutil.copy(self.groups, "/home/neil/git/lama/stats/rscripts/test_data_for_R_LM/groups.csv")
            # fit the data to a linear model and extrat the t-statistics
            # forumla is s string (eg. 'genotype' or 'genotype+sex')
            cmd = ['Rscript',
                   self.rscript,
                   voxel_file,
                   self.groups,
                   line_level_pval_out_file,
                   line_level_tstat_out_file,
                   self.formula]

            try:
                subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                msg = "R linear model failed: {}".format(e.output)
                logging.error(msg)
                raise RuntimeError("R linear model failed: {}".format(e.output))

            # Read in the pvalue and tvalue results. This will contain values from the line level call as well as
            # the speciemn-level calls and needs to be split accordingly
            p_all = np.fromfile(line_level_pval_out_file, dtype=np.float64).astype(np.float32)
            t_all = np.fromfile(line_level_tstat_out_file, dtype=np.float64).astype(np.float32)

            # Convert all NANs in the pvalues to 1.0. Need to check that this is appropriate
            p_all[np.isnan(p_all)] = 1.0

            # Convert NANs to 0. We get NAN when for eg. all input values are 0
            t_all[np.isnan(t_all)] = 0.0

            # The binary outputs from R will contain the line level call with each speciemn-level call concatenated
            # after it
            p_line = p_all[:num_data_points]
            t_line = t_all[:num_data_points]

            line_level_pvals.append(p_line)
            line_level_tvals.append(t_line)

            # Get the specimen-level statistics
            for i, row in mutants_df.iterrows():
                id_ = row['volume_id']
                start = num_data_points * (i + 1)
                end = num_data_points * (i + 2)
                print(t_all.size)
                t = t_all[start:end]
                p = p_all[start:end]
                specimen_tstats[id_].append(t)
                specimen_pvals[id_].append(p)


        line_pvals_array = np.hstack(line_level_pvals)
        line_tvals_array = np.hstack(line_level_tvals)

        try:
            os.remove(line_level_pval_out_file)
        except OSError:
            logging.info('tried to remove temporary file {}, but could not find it'.format(line_level_pval_out_file))
        try:
            os.remove(line_level_tstat_out_file)
        except OSError:
            logging.info('tried to remove temporary file {}, but could not find it'.format(line_level_tstat_out_file))


        pval_file = join(self.outdir, 'tempPvals.bin')

        line_pvals_array.tofile(pval_file)
        self.tstats = line_tvals_array
        qval_outfile = join(self.outdir, 'qvals.bin')
        pvals_distribution_image_file = join(self.outdir, PVAL_DIST_IMG_FILE)

        cmdFDR = ['Rscript',
               self.rscriptFDR,
               pval_file,
               str(line_pvals_array.shape[0]),
               qval_outfile,
               pvals_distribution_image_file
               ]

        try:
            subprocess.check_output(cmdFDR)
        except subprocess.CalledProcessError as e:
            logging.warn("R FDR calculation failed: {}".format(e.message))
            raise

        self.qvals = np.fromfile(qval_outfile, dtype=np.float64).astype(np.float32)
        self.pvals = line_pvals_array.ravel()

        self.specimen_pvals = specimen_pvals
        self.specimen_tstats = specimen_tstats


class LinearModelR(StatsTestR):
    def __init__(self, *args):
        super(LinearModelR, self).__init__(*args)
        self.rscript = join(os.path.dirname(os.path.realpath(__file__)), LINEAR_MODEL_SCRIPT)
        self.rscriptFDR = join(os.path.dirname(os.path.realpath(__file__)), PADJUST_SCRIPT)
        self.STATS_NAME = 'LinearModelR'


class CircularStatsTest(StatsTestR):
    def __init__(self, *args):
        super(CircularStatsTest, self).__init__(*args)
        self.rscript = join(os.path.dirname(os.path.realpath(__file__)), CIRCULAR_SCRIPT)
        self.STATS_NAME = 'CircularStats'
        self.MIN_DEF_MAGNITUDE = 10
        # Todo: one doing N1, can't just use Z-score as we have angles

    def run(self):
        axis = 0

        # get indices in mutants where deformation is less than a certain magnitude
        mut_magnitudes = np.linalg.norm(self.wt_data, axis=0)

        wt_bar = circmean(self.wt_data, axis=axis)
        mut_bar = circmean(self.mut_data, axis=axis)

        # Find the mutant mean that gives us the shortest distance from the WT mean
        mut_mean1 = abs(mut_bar - wt_bar)
        mut_mean2 = abs(mut_bar - (-wt_bar + 360))

        both_means = np.vstack((mut_mean1, mut_mean2))
        mut_min_mean = np.amin(both_means, axis=0)

        wt_var = circvar(self.wt_data, axis=axis)
        mut_var = circvar(self.mut_data, axis=axis)
        wt_n = len(self.wt_data)
        mut_n = len(self.mut_data)

        pvals, tstats = welch_ttest(wt_bar, wt_var, wt_n, mut_min_mean, mut_var, mut_n)
        fdr = self.fdr_class(pvals)
        qvals = fdr.get_qvalues()
        qvals[mut_magnitudes < self.MIN_DEF_MAGNITUDE] = 1.0
        self.qvals = qvals
        tstats[mut_magnitudes < self.MIN_DEF_MAGNITUDE] = 0.0
        self.tstats = tstats


class TTest(AbstractStatisticalTest):
    """
    Compare all the mutants against all the wild type. Generate a stats overlay

    TODO: Change how it calls BH as BH no longer takes a mask
    When working with
    """
    def __init__(self, *args):
        super(TTest, self).__init__(*args)
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

        return
        tstats, pvals = self.runttest(wt_chunks, mut_chunks)
        pval_chunk[np.isnan(pval_chunk)] = 1.0
        pval_chunk = pval_chunk.astype(np.float32)
        tstats.extend(tstats_chunk)
        pvals.extend(pval_chunk)

        pvals = np.array(pvals)
        tstats = np.array(tstats)

        fdr = self.fdr_class(pvals)


        qvalues = fdr.get_qvalues()
        gc.collect()

        self.filtered_tscores = self._result_cutoff_filter(tstats, qvalues) # modifies tsats in-place

        # Remove infinite values
        self.filtered_tscores[self.filtered_tscores > MINMAX_TSCORE] = MINMAX_TSCORE
        self.filtered_tscores[self.filtered_tscores < -MINMAX_TSCORE] = - MINMAX_TSCORE

    #@profile
    def runttest(self, wt, mut):

        return ttest_ind(mut, wt)

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
    def get_qvalues(self):
        """
        Mask ndarray of booleans. True == masked
        """

        # Write out pvalues to temporary file for use in R
        pvals = self.pvalues
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
        return pvals_resorted


    def ecdf(self, x):
        '''no frills empirical cdf used in fdrcorrection
        '''
        nobs = len(x)
        return np.arange(1,nobs+1)/float(nobs)


# class BenjaminiHochbergR(AbstractFalseDiscoveryCorrection):
#     def __init__(self, *args):
#         super(BenjaminiHochbergR, self).__init__(*args)
#
#     def get_qvalues(self, mask):
#         """
#         Mask ndarray of booleans. True == masked
#         """
#         print 'Doing r calculation'
#         self.pvalues[mask == False] = robj.NA_Real
#         qvals = np.array(rstats.p_adjust(FloatVector(self.pvalues), method='BH'))
#         qvals[np.isnan(qvals)] = 1
#         qvals[np.isneginf(qvals)] = 1
#         qvals[np.isinf(qvals)] = 1
#         return qvals


class Zmap(object):
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

    def process_mutant(self, mut_data, fdr=True):
        """
        Get the pixel-wise z-score of a mutant

        Parameters
        ----------
        mut_data: numpy ndarray
            1D masked array
        fdr: bool
            whether to perform fdr correction on the Z scores

        Returns
        -------
        1D np.ndarray of zscore values

        """

        z_scores = scipystats.zmap(mut_data, self.wt_data)

        # Scale inf values
        z_scores[z_scores > MINMAX_TSCORE] = MINMAX_TSCORE
        z_scores[z_scores < -MINMAX_TSCORE] = - MINMAX_TSCORE

        # Remove nans
        z_scores[np.isnan(z_scores)] = 0

        # Get probability from z-score then corret with B&H
        p = norm.sf(abs(z_scores))
        if fdr:
            logging.info("Performing BH FDR correction on Z-score pvalues")
            rejected_null = multitest.multipletests(p, method='fdr_bh')[0]  # [0] reject null @ 0.05, [1] are the q values
            z_scores[~rejected_null] = 0
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


def numpy_to_dat_1d(mat, outfile):

    # create a binary file
    binfile = file(outfile, 'wb')
    # and write out two integers with the row and column dimension

    header = struct.pack('2I', mat.shape[0])
    binfile.write(header)
    # then loop over columns and write each
    for i in range(mat.shape[1]):
        data = struct.pack('%id' % mat.shape[0], *mat[:, i])
        binfile.write(data)

    binfile.close()

