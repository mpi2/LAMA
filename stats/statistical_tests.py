import pandas as pd
import scipy.stats.stats as scipystats
from scipy.stats import norm
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
from lib import addict

MINMAX_TSCORE = 50  # If we get very large tstats or in/-inf this is our new max/min
PADJUST_SCRIPT = 'rscripts/r_padjust.R'
LINEAR_MODEL_SCRIPT = 'rscripts/lmFast.R'
CIRCULAR_SCRIPT = 'circular.R'
VOLUME_METADATA_NAME = 'volume_metadata.csv'
DATA_FILE_FOR_R_LM = 'tmp_data_for_lm'
PVAL_R_OUTFILE = 'tmp_pvals_out.dat'
TVAL_R_OUTFILE = 'tmp_tvals_out.dat'
GROUPS_FILE_FOR_LM = 'groups.csv'
STATS_FILE_SUFFIX = '_stats_'
PVAL_DIST_IMG_FILE = 'pval_distribution.png'
R_CHUNK_SIZE = 5000000


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
        self.line_tstats = None
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

        # Make temp files for storing the output of the line levl LM results
        line_level_pval_out_file = tempfile.NamedTemporaryFile().name
        line_level_tstat_out_file = tempfile.NamedTemporaryFile().name

        data = np.vstack((self.wt_data, self.mut_data))

        num_pixels = data.shape[1]

        num_chunks = num_pixels / R_CHUNK_SIZE if num_pixels > R_CHUNK_SIZE else 1

        logging.info('using formula {}'.format(self.formula))
        print('num chunks', num_chunks)

        # Loop over the data in chunks
        # np.array_split provides a split view on the array so does not increase memory
        # The result will be a bunch of arrays split across the second dimension
        chunked_data = np.array_split(data, num_chunks, axis=1)

        data_test_file = '/home/neil/Desktop/t/lama_test_df.csv'

        # These contain the chunked stats results
        line_level_pvals = []
        line_level_tvals = []

        # These contain the specimen-level statstics
        specimen_pvals = defaultdict(list)
        specimen_tstats = defaultdict(list)

        # Load in the groups file so we can get the specimen names
        groups_df = pd.read_csv(self.groups)
        mutants_df = groups_df[groups_df.genotype == 'mutant']

        i = 0  # enumerate!
        voxel_file = tempfile.NamedTemporaryFile().name

        for data_chunk in chunked_data:

            current_chink_size = data_chunk.shape[1]  # Not all chunks wil be same size
            print('chunk: {}'.format(i))
            i += 1

            numpy_to_dat(data_chunk, voxel_file)

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

            # The first chunk of data will be from the line-level call
            p_line = p_all[:current_chink_size]
            t_line = t_all[:current_chink_size]

            line_level_pvals.append(p_line)
            line_level_tvals.append(t_line)

            # Get the specimen-level statistics
            for r, (_, row) in enumerate(mutants_df.iterrows()):
                id_ = row['volume_id']
                start = current_chink_size * (r + 1)
                end = current_chink_size * (r + 2)
                t = t_all[start:end]
                p = p_all[start:end]
                specimen_tstats[id_].append(t)
                specimen_pvals[id_].append(p)

        line_pvals_array = np.hstack(line_level_pvals)
        line_tvals_array = np.hstack(line_level_tvals)

        self.line_qvals = self.do_fdr(line_pvals_array)

        try:
            os.remove(line_level_pval_out_file)
        except OSError:
            logging.info('tried to remove temporary file {}, but could not find it'.format(line_level_pval_out_file))
        try:
            os.remove(line_level_tstat_out_file)
        except OSError:
            logging.info('tried to remove temporary file {}, but could not find it'.format(line_level_tstat_out_file))

        self.line_tstats = line_tvals_array

        self.line_pvals = line_pvals_array.ravel()

        # Join up the results chunks for the specimen-level analysis. Do FDR correction on the pvalues
        self.specimen_results = addict.Dict()

        for id_, pvals in list(specimen_pvals.items()):
            p_ = np.array(pvals).ravel()
            q_ = self.do_fdr(p_)
            t_ = np.array(specimen_tstats[id_]).ravel()
            self.specimen_results[id_]['histogram'] = np.histogram(p_, bins=100)[0]
            self.specimen_results[id_]['q'] = q_
            self.specimen_results[id_]['t'] = t_
            self.specimen_results[id_]['p'] = p_

        # plt.xlim(0, 1)
        # plt.savefig(join('/home/neil/Desktop/t/t/aoefua.png'))

        # self.specimen_qvals = {id_: self.do_fdr()) for id_, pvals in  specimen_pvals.items()}
        # self.specimen_tstats = {id_: np.hstack(tstats) for id_, tstats in specimen_tstats.items()}

    def do_fdr(self, pvals):
        """
        Use R for FDR correction
        Parameters
        ----------
        pvals: numpy.ndarray
            The pvalues to be corrected
        Returns
        -------
        numpyndarray
            The corrected qvalues
        """
        pval_file = join(self.outdir, 'tempPvals.bin')
        try:
            pvals.tofile(pval_file)
        except IOError:
            pass

        qval_outfile = join(self.outdir, 'qvals.bin')
        # pvals_distribution_image_file = join(self.outdir, PVAL_DIST_IMG_FILE)

        cmdFDR = ['Rscript',
                  self.rscriptFDR,
                  pval_file,
                  str(pvals.shape[0]),
                  qval_outfile
                  ]

        try:
            subprocess.check_output(cmdFDR)
        except subprocess.CalledProcessError as e:
            logging.warn("R FDR calculation failed: {}".format(e.message))
            raise

        return np.fromfile(qval_outfile, dtype=np.float64).astype(np.float32)


class LinearModelR(StatsTestR):
    def __init__(self, *args):
        super(LinearModelR, self).__init__(*args)
        self.rscript = join(os.path.dirname(os.path.realpath(__file__)), LINEAR_MODEL_SCRIPT)
        self.rscriptFDR = join(os.path.dirname(os.path.realpath(__file__)), PADJUST_SCRIPT)
        self.STATS_NAME = 'LinearModelR'


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
    binfile = open(outfile, 'wb')
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

