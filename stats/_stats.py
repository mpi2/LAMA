import numpy as np
import scipy.stats.stats as scipystats
from scipy.special import stdtr
from scipy.stats import circmean, circvar, circstd, ttest_ind, norm
import gc
from os.path import join
import os.path
import csv
from collections import defaultdict
import subprocess
import sys
import struct
import logging
import tempfile
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




class LinearModelPython(AbstractStatisticalTest):
    def __init__(self, *args):
        super(LinearModelPython, self).__init__(*args)
        self.fdr_class = BenjaminiHochberg

        # The output of the test
        self.tstats = None
        self.qvals = None
        self.fdr_tstats = None

    def run(self):
        data = np.vstack((self.wt_data, self.mut_data))

    @staticmethod
    def cov(x, y):
        X = np.column_stack([x, y]).astype(np.float32)
        X -= X.mean(axis=0)
        N = len(y)
        fact = N
        # by_hand = np.dot(X.T, X.conj()) / fact

        x1 = X.T
        x2 = X.conj()
        x_test = np.copy(x1)
        x_test[0][0] = 4

        cov_mat = np.matmul([x1, x1, x_test], [x2, x2, x2]) / fact

        return cov_mat

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

        pval_out_file = tempfile.NamedTemporaryFile().name
        tval_out_file = tempfile.NamedTemporaryFile().name

        data = np.vstack((self.wt_data, self.mut_data))

        num_pixels = data.shape[1]

        num_chunks = num_pixels / R_CHUNK_SIZE
        if num_pixels < R_CHUNK_SIZE:
            num_chunks = 1
        logging.info('using formula {}'.format(self.formula))
        print 'num chunks', num_chunks

        # Loop over the data in chunks
        chunked_data = np.array_split(data, num_chunks, axis=1)

        #  Yaml file for quickly loading results into VPV
        # vpv_config_file = join(stats_outdir, self.output_prefix + '_VPV.yaml')
        # vpv_config = {}

        # These contain the chunked stats results
        pvals = []
        tvals = []

        i = 0
        pixel_file = tempfile.NamedTemporaryFile().name
        for data_chucnk in chunked_data:
            print 'chunk: {}'.format(i)
            i += 1

            #stacked = np.vstack(data_chucnk)
            numpy_to_dat(data_chucnk, pixel_file)

            # fit the data to a linear model and extrat the t-statistics
            # forumla is s string (eg. 'genotype' or 'genotype+sex')
            cmd = ['Rscript',
                   self.rscript,
                   pixel_file,
                   self.groups,
                   pval_out_file,
                   tval_out_file,
                   self.formula]

            try:
                subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                msg = "R linear model failed: {}".format(e.output)
                logging.error(msg)
                raise RuntimeError("R linear model failed: {}".format(e.output))

            # Read in the pvalue and tvalue results
            p = np.fromfile(pval_out_file, dtype=np.float64).astype(np.float32)
            t = np.fromfile(tval_out_file, dtype=np.float64).astype(np.float32)

            # Convert all NANs in the pvalues to 1.0. Need to check that this is appropriate
            p[np.isnan(p)] = 1.0
            pvals.append(p)

            # Convert NANs to 0. We get NAN when for eg. all input values are 0
            t[np.isnan(t)] = 0.0

            tvals.append(t)

        pvals_array = np.hstack(pvals)


        # Remove the temp data files
        # try:
        #     os.remove(pixel_file)
        # except OSError:
        #     logging.info('tried to remove temporary file {}, but could not find it'.format(pixel_file))
        try:
            os.remove(pval_out_file)
        except OSError:
            logging.info('tried to remove temporary file {}, but could not find it'.format(pval_out_file))
        try:
            os.remove(tval_out_file)
        except OSError:
            logging.info('tried to remove temporary file {}, but could not find it'.format(tval_out_file))

        tvals_array = np.hstack(tvals)
        pval_file = join(self.outdir, 'tempPvals.bin')

        # Testing
        #pvals_array = pvals_array[(tvals_array > 0) & (tvals_array < 8)]

        pvals_array.tofile(pval_file)
        self.tstats = tvals_array
        qval_outfile = join(self.outdir, 'qvals.bin')
        pvals_distribution_image_file = join(self.outdir, PVAL_DIST_IMG_FILE)

        cmdFDR = ['Rscript',
               self.rscriptFDR,
               pval_file,
               str(pvals_array.shape[0]),
               qval_outfile,
               pvals_distribution_image_file
               ]
        logging.info("Pvals array length: {}".format(pvals_array.shape[0]))

        try:
            subprocess.check_output(cmdFDR)
        except subprocess.CalledProcessError as e:
            logging.warn("R FDR calculation failed: {}".format(e.message))
            raise

        self.qvals = np.fromfile(qval_outfile, dtype=np.float64).astype(np.float32)
        self.pvals = pvals_array.ravel()

        min_t = min(tvals_array)
        try:
            t_threshold = tvals_array[(tvals_array > 0) & (self.qvals <= 0.05)].min()
        except ValueError:
            t_threshold = 'No t-statistics below fdr threshold'
        max_t = max(tvals_array)
        min_p = min(pvals_array)
        min_q = min(self.qvals)
        logging.info("\n\nMinimum T score: {}\nMaximum T score: {}\nT threshold at FDR 0.05: {}\nMinimum p-value: {}\nMinimum q-value: {}".format(
            min_t, max_t, t_threshold, min_p, min_q
        ))


        # fdr = self.fdr_class(pvals_array)
        # self.qvals = fdr.get_qvalues()


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

        z_scores = scipystats.zmap(mut_data, self.wt_data)

        # Scale inf values
        z_scores[z_scores > MINMAX_TSCORE] = MINMAX_TSCORE
        z_scores[z_scores < -MINMAX_TSCORE] = - MINMAX_TSCORE

        # Remove nans
        z_scores[np.isnan(z_scores)] = 0

        # Do fdr correction
        p = norm.sf(abs(z_scores))
        logging.info("Performing BH FDR correction on Z-score pvalues")
        corrected_significant = multitest.multipletests(p, method='fdr_bh')[0]  # [0] bool <0.05, [1] is the q values

        #filter z-scores where not significant
        z_scores[~corrected_significant] = 0
        return z_scores

class OneAgainstManytestAngular(OneAgainstManytest):
    def __init__(self, *args):
        super(OneAgainstManytestAngular, self).__init__(*args)

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
        wt_data = np.array(self.wt_data)
        angular_std = circstd(wt_data, axis=0)
        angular_mean = circmean(wt_data, axis=0)

        mut_var = mut_data - angular_mean
        mut_var[mut_var > 180] = - (360 - mut_var[mut_var > 180]) - angular_mean[mut_var > 180]

        angular_z = mut_var / angular_std

        # Filter out any values below x standard Deviations
        angular_z[np.absolute(angular_z) < self.zscore_cutoff] = 0

        # Scale inf values
        angular_z[angular_z > MINMAX_TSCORE] = MINMAX_TSCORE
        angular_z[angular_z < -MINMAX_TSCORE] = - MINMAX_TSCORE

        # Remove nans
        angular_z[np.isnan(angular_z)] = 0

        return angular_z


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

# def watson_williams(*args, **kwargs):
#     """
#     Taken from: https://github.com/circstat/pycircstat/blob/master/pycircstat/tests.py
#
#     Parametric Watson-Williams multi-sample test for equal means. Can be
#     used as a one-way ANOVA test for circular data.
#     H0: the s populations have equal means
#     HA: the s populations have unequal means
#     Note:
#     Use with binned data is only advisable if binning is finer than 10 deg.
#     In this case, alpha is assumed to correspond
#     to bin centers.
#     The Watson-Williams two-sample test assumes underlying von-Mises
#     distributrions. All groups are assumed to have a common concentration
#     parameter k.
#     :param args: number of arrays containing the data; angles in radians
#     :param w:    list the same size as the number of args containing the number of
#                  incidences for each arg. Must be passed as keyword argument.
#     :param axis: the test will be performed along this axis. Must be passed as keyword
#                  argument.
#     :return pval, table: p-value and pandas dataframe containing the ANOVA table
#     """
#
#     axis = kwargs.get('axis', None)
#     w = kwargs.get('w', None)
#
#     # argument checking
#     if w is not None:
#         assert len(w) == len(
#             args), "w must have the same length as number of arrays"
#         for i, (ww, alpha) in enumerate(zip(w, args)):
#             assert ww.shape == alpha.shape, "w[%i] and argument %i must have same shape" % (
#                 i, i)
#     else:
#         w = [np.ones_like(a) for a in args]
#
#     if axis is None:
#         alpha = list(map(np.ravel, args))
#         w = list(map(np.ravel, w))
#     else:
#         alpha = args
#
#     k = len(args)
#
#     # np.asarray(list())
#     ni = list(map(lambda x: np.sum(x, axis=axis), w))
#     ri = np.asarray([descriptive.resultant_vector_length(
#         a, ww, axis=axis) for a, ww in zip(alpha, w)])
#
#     r = descriptive.resultant_vector_length(
#         np.concatenate(
#             alpha, axis=axis), np.concatenate(
#             w, axis=axis), axis=axis)
#     # this must not be the numpy sum since the arrays are to be summed
#     n = sum(ni)
#
#     rw = sum([rii * nii / n for rii, nii in zip(ri, ni)])
#     kk = kappa(rw[None, ...], axis=0)
#
#     beta = 1 + 3. / (8 * kk)
#     A = sum([rii * nii for rii, nii in zip(ri, ni)]) - r * n
#     B = n - sum([rii * nii for rii, nii in zip(ri, ni)])
#
#     F = (beta * (n - k) * A / (k - 1) / B).squeeze()
#     pval = stats.f.sf(F, k - 1, n - k).squeeze()
#
#     if np.any((n >= 11) & (rw < .45)):
#         logging.warn(
#             'Test not applicable. Average resultant vector length < 0.45.')
#     elif np.any((n < 11) & (n >= 7) & (rw < .5)):
#         logging.warn(
#             'Test not applicable. Average number of samples per population 6 < x < 11 '
#             'and average resultant vector length < 0.5.')
#     elif np.any((n >= 5) & (n < 7) & (rw < .55)):
#         logging.warn(
#             'Test not applicable. Average number of samples per population 4 < x < 7 and '
#             'average resultant vector length < 0.55.')
#     elif np.any(n < 5):
#         logging.warn(
#             'Test not applicable. Average number of samples per population < 5.')
#
#     if np.prod(pval.shape) > 1:
#         T = np.zeros_like(pval, dtype=object)
#         for idx, p in np.ndenumerate(pval):
#             T[idx] = pd.DataFrame({'Source': ['Columns', 'Residual', 'Total'],
#                                    'df': [k - 1, n[idx] - k, n[idx] - 1],
#                                    'SS': [A[idx], B[idx], A[idx] + B[idx]],
#                                    'MS': [A[idx] / (k - 1), B[idx] / (n[idx] - k), np.NaN],
#                                    'F': [F[idx], np.NaN, np.NaN],
#                                    'p-value': [p, np.NaN, np.NaN]}).set_index('Source')
#
#     else:
#         T = pd.DataFrame({'Source': ['Columns', 'Residual', 'Total'],
#                           'df': [k - 1, n - k, n - 1],
#                           'SS': [A, B, A + B],
#                           'MS': [A / (k - 1), B / (n - k), np.NaN],
#                           'F': [F, np.NaN, np.NaN],
#                           'p-value': [pval, np.NaN, np.NaN]}).set_index('Source')
#
#     return pval, T
#
#
#

def welch_ttest(abar, avar, na, bbar, bvar, nb):
    adof = na - 1
    bdof = nb - 1
    tf = (abar - bbar) / np.sqrt(avar/na + bvar/nb)
    dof = (avar/na + bvar/nb)**2 / (avar**2/(na**2*adof) + bvar**2/(nb**2*bdof))
    pf = 2*stdtr(dof, -np.abs(tf))
    return pf, tf