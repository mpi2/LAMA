import numpy as np
import scipy.stats.stats as stats
import gc
from os.path import join
import os.path
import csv
from collections import defaultdict
import subprocess
import sys
import struct
import logging

import SimpleITK as sitk

sys.path.insert(0, join(os.path.dirname(__file__), '..'))

MINMAX_TSCORE = 50 # If we get very large tstats or in/-inf this is our new max/min
# PADJUST_SCRIPT = 'r_padjust.R'
LINEAR_MODEL_SCIPT = 'lmFast.R'
FDR_SCRPT = 'r_padjust.R'
VOLUME_METADATA_NAME = 'volume_metadata.csv'
DATA_FILE_FOR_R_LM = 'tmp_data_for_lm'
PVAL_R_OUTFILE = 'tmp_pvals_out.dat'
TVAL_R_OUTFILE = 'tmp_tvals_out.dat'
GROUPS_FILE_FOR_LM = 'groups.csv'
STATS_FILE_SUFFIX = '_stats_'
FDR_CUTOFF = 0.05


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

    def get_volume_metadata(self):
        """
        Get the metada for the volumes, such as sex and (in the future) scan date
        Not currently used. Superceded by method in _reg_stats_new.py
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

    def write_result(self, result_array, outpath):
        """
        """
         # Create a full size output array
        size = np.prod(self.shape)
        full_output = np.zeros(size)

        # Insert the result p and t vals back into full size array
        full_output[self.mask != False] = result_array

        reshaped_results = full_output.reshape(self.shape)
        result_img = sitk.GetImageFromArray(reshaped_results)
        sitk.WriteImage(result_img, outpath, True)


class LinearModelR(AbstractStatisticalTest):
    def __init__(self, *args):
        super(LinearModelR, self).__init__(*args)
        self.stats_method_object = None  # ?
        self.fdr_class = BenjaminiHochberg

        self.STATS_NAME = 'LinearModelR'

        self.tstats = None
        self.qvals = None
        self.fdr_tstats = None

    def set_formula(self, formula):
        self.formula = formula

    def set_groups(self, groups):
        self.groups = groups

    def run(self):

        if not self.groups:
            # We need groups file for linera model
            logging.warn('linear model failed. We need groups file')
            return

        # np.array_split provides a split view on the array so does not increase memory
        # The result will be a bunch of arrays split across the second dimension

        linear_model_script = join(os.path.dirname(os.path.realpath(__file__)), LINEAR_MODEL_SCIPT)

        pval_out_file = join(self.outdir, PVAL_R_OUTFILE)
        tval_out_file = join(self.outdir, TVAL_R_OUTFILE)

        # TODO: this needs changing as it takes too much memory
        data = np.vstack((self.wt_data, self.mut_data))

        num_pixels = data.shape[1]
        chunk_size = 200000
        num_chunks = num_pixels / chunk_size
        if num_pixels < 200000:
            num_chunks = 1
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
        for data_chucnk in chunked_data:
            i += 1
            pixel_file = join(self.outdir, DATA_FILE_FOR_R_LM)
            numpy_to_dat(np.vstack(data_chucnk), pixel_file)

            # fit the data to a linear model and extrat the tvalue
            cmd = ['Rscript',
                   linear_model_script,
                   pixel_file,
                   self.groups,
                   pval_out_file,
                   tval_out_file,
                   self.formula]

            try:
                subprocess.check_output(cmd)
            except subprocess.CalledProcessError as e:
                logging.warn("R linear model failed: {}".format(e))
                raise

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
        try:
            os.remove(pixel_file)
        except OSError:
            logging.info('tried to remove temporary file {}, but could not find it'.format(pixel_file))
        try:
            os.remove(pval_out_file)
        except OSError:
            logging.info('tried to remove temporary file {}, but could not find it'.format(pval_out_file))
        try:
            os.remove(tval_out_file)
        except OSError:
            logging.info('tried to remove temporary file {}, but could not find it'.format(tval_out_file))

        tvals_array = np.hstack(tvals)

        self.tstats = tvals_array
        fdr = self.fdr_class(pvals_array)
        self.qvals = fdr.get_qvalues()

class TTest(AbstractStatisticalTest):
    """
    Compare all the mutants against all the wild type. Generate a stats overlay

    TODO: Change how it calls BH as BH no longer takes a mask
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
        qvalues = fdr.get_qvalues()
        gc.collect()

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

        z_scores = stats.zmap(mut_data, self.wt_data)

        # Filter out any values below x standard Deviations
        z_scores[np.absolute(z_scores) < self.zscore_cutoff] = 0

        # Scale inf values
        z_scores[z_scores > MINMAX_TSCORE] = MINMAX_TSCORE
        z_scores[z_scores < -MINMAX_TSCORE] = - MINMAX_TSCORE

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
