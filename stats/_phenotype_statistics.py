#!/usr/bin/env python

"""

"""

import os
import sys
from os.path import join, basename, split, isdir, isfile

# Hack. Relative package imports won't work if this module is run as __main__
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
from lib import addict
import common
import SimpleITK as sitk
from elastix.invert import InvertSingleVol, InvertStats
from statistical_tests import OneAgainstManytest, OneAgainstManytestAngular
from data_getters import DeformationDataGetter, IntensityDataGetter, JacobianDataGetter, AngularDataGetter
import numpy as np
import gc
from statistical_tests import LinearModelR, CircularStatsTest
import logging
import shutil
import tsne
from automated_annotation import Annotator
import csv
import pandas as pd
from scipy.stats import zmap

STATS_FILE_SUFFIX = '_stats_'
CALC_VOL_R_FILE = 'rscripts/calc_organ_vols.R'
CLUSTER_PLOT_NAME = 'n1_clustering.png'
MINMAX_TSCORE = 50
FDR_CUTOFF = 0.05


class AbstractPhenotypeStatistics(object):
    """
    The base class for the statistics generators
    """
    def __init__(self, out_dir, analysis_prefix, main_config, analysis_config):
        """
        Parameters
        ----------
        mask_array: numpy ndarray
            1D mask array
        groups: dict
            specifies which groups the data volumes belong to (for linear model etc.)
        label_map: ndarray
            labels for calculating organ volumes
        label_names: Dict
            {0: 'organ name1', 1: 'organ name2' ....}

        """
        self.blur_fwhm = main_config.blur_fwhm
        self.root_dir = main_config.root_dir
        self.normalisation_roi = main_config.normalisation_roi
        self.subsampled_mask = None # Not been using that so deprecate it
        self.n1 = main_config.n1
        self.label_map = main_config.label_map
        self.label_names = main_config.label_names
        self.project_name = main_config.project_name
        self.out_dir = out_dir
        common.mkdir_if_not_exists(self.out_dir)
        self.mask = main_config.mask_array_flat  # this is a flat binary array
        self.formulas = main_config.formulas
        self.wt_file_list = main_config.wt_file_list
        self.mut_file_list = main_config.mut_file_list
        self.voxel_size = main_config.voxel_size
        self.n1_out_dir = join(self.out_dir, 'n1')
        self.filtered_stats_path = None
        self.stats_out_dir = None
        self.n1_tester = OneAgainstManytest
        self.analysis_config = analysis_config

        # Obtained from the datagetter
        self.shape = None

        self.n1_stats_output = []  # Paths to the n1 anlaysis output. Use din inverting stats volumes
        self.groups = main_config.groups


    def _set_data(self):
        """
        Set the wt and mut data.
        """

        vol_order = self.get_volume_order()
        self.dg = self.data_getter(self.wt_file_list, self.mut_file_list, self.mask, vol_order, self.voxel_size,
                                    self.subsampled_mask, None, self.blur_fwhm)  # None is subsample which is deprecated

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
                    if first:  # Skip header
                        first = False
                        continue
                    else:
                        order.append(row[0])
            return order
        else:
            return None

    def run(self, stats_object, analysis_prefix):

        self.analysis_prefix = analysis_prefix
        try:
            self._set_data()
        except IOError as e:
            print('error getting data for {}: {}'.format(self.analysis_prefix, e))
            return False
        normalisation_dir = join(self.out_dir, 'normalised_images')
        self.dg.set_normalisation_roi(self.normalisation_roi, normalisation_dir)  # only used for ItensityStats
        self.dg.set_data()

        logging.info('using wt_paths n={}\n--------------\n{}\n\n'.format(
            len(self.dg.wt_paths), '\n'.join([x for x in self.dg.wt_paths])))

        logging.info('using mut_paths n={}\n--------------\n{}\n\n'.format(
            len(self.dg.mut_paths), '\n'.join([x for x in self.dg.mut_paths])))

        self.shape = self.dg.shape
        self._many_against_many(stats_object)
        if self.n1:
            self._one_against_many()
        del self.dg
        gc.collect()

    def _one_against_many(self):
        """
        Compare each mutant seperatley against all wildtypes
        """
        n1 = self.n1_tester(self.dg.masked_wt_data)
        common.mkdir_if_not_exists(self.n1_out_dir)

        self.n1_prefix = self.analysis_prefix + STATS_FILE_SUFFIX

        for path, mut_data in zip(self.dg.mut_paths, self.dg.masked_mut_data):
            result = n1.process_mutant(mut_data)
            reshaped_data = np.zeros(np.prod(self.shape))
            reshaped_data[self.mask != False] = result
            reshaped_data = reshaped_data.reshape(self.shape)
            out_path = join(self.n1_out_dir, self.n1_prefix + os.path.basename(path))
            self.n1_stats_output.append(out_path)
            common.write_array(reshaped_data, out_path)
        del n1
        # Do some clustering on the Zscore results in order to identify potential partial penetrence
        tsne_plot_path = join(self.out_dir, CLUSTER_PLOT_NAME)
        try:
            tsne_labels = tsne.cluster(self.n1_out_dir, tsne_plot_path)
        except (ValueError, AssertionError): # sometimes fails. Think it might be when images are identical during tests
            logging.warning('t-sne clustering failed')
        else:
            labels_str = "\n***clustering plot labels***\n"
            for num, name in tsne_labels.iteritems():
                labels_str += "{}: {}\n".format(num, name)
            logging.info(labels_str)
        gc.collect()

    def _many_against_many(self, stats_object):
        """
        Compare all mutants against all wild types
        """

        for formula in self.formulas[:1]:  # Just do one formula for now as it may break
            so = stats_object(self.dg.masked_wt_data, self.dg.masked_mut_data, self.shape, self.out_dir)

            if type(so) in (LinearModelR, CircularStatsTest):
                logging.info(common.git_log())
                so.set_formula(formula)
                so.set_groups(self.groups)
                so.run()
                qvals = so.qvals
                tstats = so.tstats
                pvals = so.pvals
                unmasked_tstats = self.rebuid_masked_output(tstats, self.mask, self.mask.shape).reshape(self.shape)
                unmasked_qvals = self.rebuid_masked_output(qvals, self.mask, self.mask.shape).reshape(self.shape)
                unmasked_pvals = self.rebuid_masked_output(pvals, self.mask, self.mask.shape).reshape(self.shape)
                filtered_tsats = self.write_results(unmasked_qvals, unmasked_tstats,  unmasked_pvals, so.STATS_NAME, formula)
                t_threshold_file = join(self.out_dir, 'Qvals-{}.csv'.format(self.type))
                write_threshold_file(unmasked_qvals, unmasked_tstats, t_threshold_file)

                del so
                gc.collect()
            else:
                so.run()
                qvals = so.qvals
                tstats = so.tstats
                fdr_tsats = so.fdr_tstats
                filtered_tsats = self.write_results(qvals, tstats, fdr_tsats, self.mask)
                del so
            if self.label_map is not None and self.label_names:
                logging.info("Doing auto annotation")
                ann_outpath = join(self.out_dir, 'annotation.csv')
                ann = Annotator(self.label_map, self.label_names, filtered_tsats, ann_outpath)
                ann.annotate()
            else:
                logging.info("Skipping auto annotation as there was either no labelmap or list of label names")

    def rebuid_masked_output(self, array, mask, shape):
        """
        The results from the stats objects have masked regions removed. Add the result back into a full-sized image
        Override this method for subsampled analysis e.g. GLCM
        """
        array[array > MINMAX_TSCORE] = MINMAX_TSCORE
        array[array < -MINMAX_TSCORE] = - MINMAX_TSCORE
        full_output = np.zeros(shape)
        full_output[mask != False] = array
        return full_output.reshape(shape)

    def write_results(self, qvals, tstats, pvals, stats_name, formula=None):
        # Write out the unfiltered t values and p values

        stats_prefix = self.project_name + '_' + self.analysis_prefix
        if formula:
            stats_prefix += '_' + formula
        stats_outdir = join(self.out_dir, stats_name)
        common.mkdir_if_not_exists(stats_outdir)
        unfilt_tq_values_path = join(stats_outdir,  stats_prefix + '_' + stats_name + '_t_q_stats')

        np.savez_compressed(unfilt_tq_values_path,
                            tvals=[tstats],
                            qvals=[qvals]
                            )

        self.stats_out_dir = stats_outdir
        outpath = join(stats_outdir, stats_prefix + '_' + stats_name + '_' + formula + '_FDR_' + str(0.5) + '_stats_.nrrd')
        outpath_unfiltered_tstats = join(stats_outdir, stats_prefix + '_' + stats_name + '_Tstats_' + formula + '_stats_.nrrd')
        outpath_unfiltered_pvals = join(stats_outdir, stats_prefix + '_' + stats_name + '_pvals_' + formula + '_stats_.nrrd')

        self.filtered_stats_path = outpath
        common.write_array(tstats, outpath_unfiltered_tstats)

        common.write_array(pvals, outpath_unfiltered_pvals)

        # Write filtered tstats overlay. Done here so we don't have filtered and unfiltered tstats in memory
        # at the same time
        try:
            filtered_tsats = self._result_cutoff_filter(tstats, qvals)
        except ValueError:
            print "Tstats and qvalues are not equal size"
        else:
            common.write_array(filtered_tsats, outpath)
        gc.collect()
        return filtered_tsats # The fdr-corrected stats


    def rebuid_subsamlped_output(self, array, shape, chunk_size):
        """

        Parameters
        ----------
        array: numpy.ndarray
            the subsampled array to rebuild
        shape: tuple
            the shape of the final result
        chunk_size: int
            the original subsampling factor

        Returns
        -------un
        np.ndarray
            rebuilt array of the same size of the original inputs data

        """
        out_array = np.zeros(self.shape)
        i = 0
        for z in range(0, shape[0] - chunk_size, chunk_size):
            for y in range(0, shape[1] - chunk_size, chunk_size):
                for x in range(0, shape[2] - chunk_size, chunk_size):
                    out_array[z: z + chunk_size, y: y + chunk_size, x: x + chunk_size] = array[i]
                    i += 1

        return out_array

    @staticmethod
    def _result_cutoff_filter(t, q):
        """
        Convert to numpy arrays and set to zero any tscore that has a corresponding pvalue > 0.05

        Parameters
        ----------

        """
        if len(t) != len(q):
            raise ValueError
        else:
            mask = q > FDR_CUTOFF
            t[mask] = 0

        return t

    def invert(self, invert_config_path):
        """
        Invert the stats back onto the rigidly aligned volumes

        Parameters
        ----------
        invert_order: dict
            Contains inversion order information
        """
        # TODO.

        # Invert the n1 stats
        n1_invert_out_dir = join(self.n1_out_dir, 'inverted')
        common.mkdir_if_not_exists(n1_invert_out_dir)
        for stats_vol_path in self.n1_stats_output:
            n1_inverted_out = join(n1_invert_out_dir, basename(stats_vol_path))
            inv = InvertSingleVol(invert_config_path, stats_vol_path, n1_inverted_out)
            inv.run(prefix=self.n1_prefix)

        # Invert the Linear model/ttest stats
        if not self.filtered_stats_path:
            return

        # inverted_stats
        stats_invert_dir = join(self.stats_out_dir, 'inverted')
        common.mkdir_if_not_exists(stats_invert_dir)
        invs = InvertStats(invert_config_path, self.filtered_stats_path, stats_invert_dir)
        invs.run()


class IntensityStats(AbstractPhenotypeStatistics):
    def __init__(self, *args, **kwargs):
        super(IntensityStats, self).__init__(*args, **kwargs)
        self.data_getter = IntensityDataGetter
        self.type = 'intensity'

class AngularStats(AbstractPhenotypeStatistics):
    def __init__(self, *args, **kwargs):
        super(AngularStats, self).__init__(*args, **kwargs)
        self.data_getter = AngularDataGetter
        self.n1_tester = OneAgainstManytestAngular
        self.type = 'circular'


class JacobianStats(AbstractPhenotypeStatistics):
    # Not used. Intensity and jacoabian analysis is the same
    def __init__(self, *args, **kwargs):
        super(JacobianStats, self).__init__(*args, **kwargs)
        self.data_getter = JacobianDataGetter
        self.type = 'jacobian'


class DeformationStats(AbstractPhenotypeStatistics):
    def __init__(self, *args, **kwargs):
        super(DeformationStats, self).__init__(*args, **kwargs)
        self.data_getter = DeformationDataGetter
        self.type = 'deformations'


class OrganVolumeStats(AbstractPhenotypeStatistics):
    """
    The volume organ data does not fit with the other classes above which all work at the pixel not label level
    """
    def __init__(self,  *args, **kwargs):
        super(OrganVolumeStats, self).__init__(*args, **kwargs)
        self.type = 'organvolume'

    def invert(self, _):
        """
        Override the parent invert. Organ volume stats works on already inverted volumes

        """
        return

    def run(self, stats_method_object, analysis_prefix):
        """
        TODO: Catch if len labels != dat length
        Parameters
        ----------
        stats_method_object
        analysis_prefix

        Returns
        -------

        Notes
        -----
        wt_vols_df and mut_vols_df are pd.DataFrames with volumes in columns (minus extension) and organ volumes in rows

        """
        def normalize_to_whole_embryo(label_vol_df, mask_vol_df):
            for col in label_vol_df:
                label_vol_df[col] = label_vol_df[col].div(mask_vol_df[col])

        if not self.label_names:
            logging.error('No label names csv path specified in stats.yaml config')
            return
        if self.label_map is None:
            logging.error('No label map image path specified in stats.yaml config')
            return
        labels = self.label_names.values()[1:]
        common.mkdir_if_not_exists(self.out_dir)

        wt_vols_df = self.get_label_vols(self.wt_file_list)
        mut_vols_df = self.get_label_vols(self.mut_file_list)

        # Get the   inverted masks if available
        wt_inveretd_mask_dir = self.analysis_config.get('wt_inverted_masks')
        mut_inveretd_mask_dir = self.analysis_config.get('mut_inverted_masks')


        if wt_inveretd_mask_dir and mut_inveretd_mask_dir:

            wt_inv_mask_path = join(self.root_dir, wt_inveretd_mask_dir)
            mut_inv_mask_path = join(self.root_dir, mut_inveretd_mask_dir)

            wt_inv_masks = common.get_file_paths(wt_inv_mask_path)
            mut_inv_masks = common.get_file_paths(mut_inv_mask_path)

            #DataFrames where the first row is the count of 1s (the whole embryo volume)
            wt_inv_masks_df = self.get_label_vols(wt_inv_masks).iloc[0]
            mut_inv_masks_df = self.get_label_vols(mut_inv_masks).iloc[0]

            normalize_to_whole_embryo(wt_vols_df, wt_inv_masks_df)
            normalize_to_whole_embryo(mut_vols_df, mut_inv_masks_df)

        # Raw organ volumes
        organ_volumes_path = join(self.out_dir, 'organ_volumes.csv')

        # create pandas DataFrames where rows=samples, columns=labels/organs
        mut = mut_vols_df.T
        wt = wt_vols_df.T

        if not mut.shape[1] == len(labels):
            msg = "The number of labels ({}) in the mutant label volumes does not match the number of labels ({}) in the lables name file".format(
                    mut.shape[1], len(labels))
            logging.error(msg)
            raise ValueError(msg)

        if not wt.shape[1] == len(labels):
            msg = "The number of labels ({}) in the wild type label volumes does not match the number of labels ({}) in the lables name file".format(
                    wt.shape[1], len(labels))
            logging.error(msg)
            raise ValueError(msg)


        muts_and_wts = pd.concat([mut, wt])
        muts_and_wts.columns = labels
        muts_and_wts.to_csv(organ_volumes_path)
        mut_vals = mut.values
        wt_vals = wt.values
        so = LinearModelR(wt_vals, mut_vals,  self.shape, self.out_dir)

        so.set_formula(self.formulas[0])
        so.set_groups(self.groups)
        so.run()
        qvals = so.qvals

        tstats = so.tstats
        pvals = so.pvals

        significant = ['yes'if x <= 0.05 else 'no' for x in qvals]
        volume_stats_path = join(self.out_dir, 'inverted_organ_volumes_LinearModel_FDR5%.csv')
        columns = ['p', 'q', 't', 'significant']
        stats_df = pd.DataFrame(index=labels, columns=columns)
        stats_df['p'] = pvals
        stats_df['q'] = qvals
        stats_df['t'] = tstats
        stats_df['significant'] = significant
        stats_df = stats_df.sort_values('q')
        stats_df.to_csv(volume_stats_path)

        # Z-scores l
        #Swith off for now. Problem is that labels is a list of dicts and needs to be a list of organ names
        # zscore_stats_path = join(self.out_dir, 'organ_volume_z_scores.csv')
        # zscores = zmap(mut_vols_df.T, wt_vols_df.T)
        # specimens = mut_vols_df.columns
        # z_df = pd.DataFrame(index=specimens, columns=labels)
        # z_df[:] = zscores
        # z_df.to_csv(zscore_stats_path)

    @staticmethod
    def get_label_vols(label_paths):
        """

        Parameters
        ----------
        label_paths: str
            paths to labelmap volumes

        Returns
        -------
        Dict: {volname:label_num: [num_voxels_1, num_voxels2...]...}
        """

        label_volumes = addict.Dict()
        for label_path in label_paths:
            # Get the name of the volume
            volname = os.path.split(split(label_path)[0])[1]
            labelmap = sitk.ReadImage(label_path)

            lsf = sitk.LabelStatisticsImageFilter()
            labelmap = sitk.Cast(labelmap, sitk.sitkUInt16)
            lsf.Execute(labelmap, labelmap)
            num_labels = lsf.GetNumberOfLabels()
            for i in range(1, num_labels + 1):
                voxel_count= lsf.GetCount(i)
                label_volumes[volname][i] = voxel_count
        return pd.DataFrame(label_volumes.to_dict())


def write_threshold_file(pvals, tvals, outpath):
    """
    Replicate the 'Qvlas-intensities/jacobians.csv' output by the TCP pipeline. Needed for gettting our data up onto
    IMPC pipeline. All we neeed is the first and last column, so just ad 'NA' for all the others

        An eample file looks like this
            "","F-statistic","tvalue-(Intercept)","tvalue-gf$genotypeKO"
            "0.01",NA,6.04551674399839,NA
            "0.05",NA,4.063447298063,NA
            "0.1",30.8843220744942,3.27694469338307,5.55736646933547
            "0.15",20.2650883331287,2.83426768232588,4.50167616928725
            "0.2",15.2004082182636,2.51876041070957,3.89877009045976
    """

    rows = ['"","F-statistic","tvalue-(Intercept)","tvalue-gf$genotypeKO"\n']
    row_template = '"{}", NA, NA, {}\n'
    for pvalue in [0.000001, 0.00001, 0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2]:
        try:
            t_thresh = np.min(tvals[np.where((pvals <= pvalue) & (tvals > 0))])
        except ValueError:  # No minimum availbale
            t_thresh = 'NA'

        row = row_template.format(str(pvalue), str(t_thresh))
        rows.append(row)
    with open(outpath, 'w') as fh:
        for r in rows:
            fh.write(r)


if __name__ == '__main__':
    # Just dev out the p/t threshold file
    import sys
    pt = sys.argv[1]  # lama npz file containg q and t values
    out = sys.argv[2]

    data = np.load(pt)
    q = data['qvals'][0].astype(np.float16)
    t = data['tvals'][0].astype(np.float16)

    write_threshold_file(q, t, out)








