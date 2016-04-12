
from os.path import join, basename
import sys
import os
# Hack. Relative package imports won't work if this module is run as __main__
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common
import SimpleITK as sitk
from invert import InvertSingleVol, InvertStats
from _stats import OneAgainstManytest, OneAgainstManytestAngular
from _data_getters import GlcmDataGetter, DeformationDataGetter, ScalarDataGetter, JacobianDataGetter, AngularDataGetter
import numpy as np
import gc
import csv
import yaml
from _stats import LinearModelR, CircularStatsTest
import logging

STATS_FILE_SUFFIX = '_stats_'
CALC_VOL_R_FILE = 'calc_organ_vols.R'
MINMAX_TSCORE = 50
FDR_CUTOFF = 0.05


class AbstractPhenotypeStatistics(object):
    """
    The base class for the statistics generators
    """
    def __init__(self, out_dir, wt_data_dir, mut_data_dir, project_name, mask_array=None, groups=None,
                 formulas=None, n1=True, voxel_size=None, wt_subset=None, mut_subset=None):
        """
        Parameters
        ----------
        mask_array: numpy ndarray
            1D mask array
        groups: dict
            specifies which groups the data volumes belong to (for linear model etc.)
        """
        self.wt_subset = wt_subset
        self.mut_subset = mut_subset
        self.n1 = n1
        self.project_name = project_name
        self.out_dir = out_dir
        common.mkdir_if_not_exists(self.out_dir)
        self.mask = mask_array  # this is a flat binary array
        self.groups = groups
        self.formulas = formulas
        self._wt_data_dir = wt_data_dir
        self._mut_data_dir = mut_data_dir
        self.voxel_size = voxel_size
        self.n1_out_dir = join(self.out_dir, 'n1')
        self.filtered_stats_path = None
        self.stats_out_dir = None
        self.n1_tester = OneAgainstManytest

        # Obtained from the datagetter
        self.shape = None

        self.n1_stats_output = []  # Paths to the n1 anlaysis output. Use din inverting stats volumes

    def _set_data(self):
        """
        Set the wt and mut data. What are the types?
        """

        vol_order = self.get_volume_order()
        self.dg = dg = self.data_getter(self._wt_data_dir, self._mut_data_dir, self.mask, vol_order, self.voxel_size,
                                        self.wt_subset, self.mut_subset)
        logging.info('using wt_paths\n--------------\n{}\n\n'.format(
            '\n'.join([basename(x) for x in self.dg.wt_paths])))

        logging.info('using mut_paths\n--------------\n{}\n\n'.format(
            '\n'.join([basename(x) for x in self.dg.mut_paths])))

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
            print 'error getting data for {}: {}'.format(self.analysis_prefix, e)
            return False
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
            outimg = sitk.GetImageFromArray(reshaped_data)
            sitk.WriteImage(outimg, out_path, True)  # Compress output
        del n1
        gc.collect()

    def _many_against_many(self, stats_object):
        """
        Compare all mutants against all wild types
        """
        so = stats_object(self.dg.masked_wt_data, self.dg.masked_mut_data, self.shape, self.out_dir)

        if type(so) in(LinearModelR, CircularStatsTest):
            for formula in self.formulas:
                so.set_formula(formula)
                so.set_groups(self.groups)
            so.run()
            qvals = so.qvals
            tstats = so.tstats
            self.write_results(qvals, tstats, so.STATS_NAME, formula)
            gc.collect()
        else:
            so.run()
            qvals = so.qvals
            tstats = so.tstats
            fdr_tsats = so.fdr_tstats
            self.write_results(qvals, tstats, fdr_tsats)
        del so

    def rebuid_output(self, array):
        """
        The results from the stats objects have masked regions removed. Add the result back into a full-sized image
        Override this method for subsampled analysis e.g. GLCM
        """
        array[array > MINMAX_TSCORE] = MINMAX_TSCORE
        array[array < -MINMAX_TSCORE] = - MINMAX_TSCORE
        full_output = np.zeros(np.prod(self.shape))
        full_output[self.mask != False] = array
        return full_output.reshape(self.shape)

    def write_results(self, qvals, tstats, stats_name, formula=None):
        # Write out the unfiltered t values and p values
        qvals = self.rebuid_output(qvals)
        tstats = self.rebuid_output(tstats)

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
        self.filtered_stats_path = outpath

        # Write filtered tstats overlay. Done here so we don't have filtered and unfiltered tstats in memory
        # at the same time
        try:
            filtered_tsats = self._result_cutoff_filter(tstats, qvals)
        except ValueError:
            print "Tstats and qvalues are not equal size"
        else:
            sitk.WriteImage(sitk.GetImageFromArray(filtered_tsats), outpath, True)
        gc.collect()

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
    def __init__(self, *args):
        super(IntensityStats, self).__init__(*args)
        self.data_getter = ScalarDataGetter


class AngularStats(AbstractPhenotypeStatistics):
    def __init__(self, *args):
        super(AngularStats, self).__init__(*args)
        self.data_getter = AngularDataGetter
        self.n1_tester = OneAgainstManytestAngular


class GlcmStats(AbstractPhenotypeStatistics):
    def __init__(self, *args):
        super(GlcmStats, self).__init__(*args)
        self.data_getter = GlcmDataGetter  # Currently just gets inertia feature with ITK default settings
        self.mask = self.create_subsampled_mask()

    def _one_against_many(self):
        """
        Not currently working
        """
        logging.info('n1 analysis not currently implemented for GLCMs')

    def _set_data(self):
        """
        Temp: Overided as we do not want shape set in this manner. Rewrite!
        """

        vol_order = self.get_volume_order()
        self.dg = self.data_getter(self._wt_data_dir, self._mut_data_dir, self.mask, vol_order)

    def get_glcm_config_values(self):
        """
        Extract glcm metadata from the glcm output folder
        """
        config_path = join(self._wt_data_dir, 'glcm.yaml')
        with open(config_path) as fh:
            config = yaml.load(fh)

        chunk_size = config['chunksize']
        original_size = config['original_shape']

        return chunk_size, original_size

    def rebuid_output(self, array):
        array[array > MINMAX_TSCORE] = MINMAX_TSCORE
        array[array < -MINMAX_TSCORE] = - MINMAX_TSCORE

        shape = self.shape # Where is  this set?
        chunk_size = self.chunk_size
        out_array = np.zeros(self.shape)
        i = 0
        for x in range(0, shape[2] - chunk_size, chunk_size):
            for y in range(0, shape[1] - chunk_size, chunk_size):
                for z in range(0, shape[0] - chunk_size, chunk_size):
                    out_array[z: z + chunk_size, y: y + chunk_size, x: x + chunk_size] = array[i]
                    i += 1

        return out_array

    def create_subsampled_mask(self):
        """
        As the glcm data is subsampled, we need a subsampled mask
        """
        chunk_size, shape = self.get_glcm_config_values()
        self.shape = shape  # This is set here as it would be the size of the subsampled glcm output
        self.chunk_size = chunk_size
        out_array = np.zeros(shape)
        i = 0
        subsampled_mask = []
        # We go x-y-z as thats how it comes out of the GLCM generator
        for x in range(0, shape[2] - chunk_size, chunk_size):
            for y in range(0, shape[1] - chunk_size, chunk_size):
                for z in range(0, shape[0] - chunk_size, chunk_size):
                    mask_region = self.mask[z: z + chunk_size, y: y + chunk_size, x: x + chunk_size]
                    if np.any(mask_region):
                        subsampled_mask.insert(i, 0)
                    else:
                        subsampled_mask.insert(i, 1)
                    i += 1

        return out_array


class JacobianStats(AbstractPhenotypeStatistics):
    # Not used. Intensity and jacoabian analysis is the same
    def __init__(self, *args):
        super(JacobianStats, self).__init__(*args)
        self.data_getter = JacobianDataGetter


class DeformationStats(AbstractPhenotypeStatistics):
    def __init__(self, *args):
        super(DeformationStats, self).__init__(*args)
        self.data_getter = DeformationDataGetter


class OrganVolumeStats(object):
    """
    The volume organ data does not fit with the other classes above
    """
    def __init__(self, outdir, wt_path, mut_path):
        self.outdir = outdir
        self.wt_path = wt_path
        self.mut_path = mut_path

    def run(self, stats_object, analysis_prefix):
        wt_data = self._get_label_vols(self.wt_path)
        mut_data = self._get_label_vols(self.mut_path)


    def _write_results(self, results, outfile):
        with open(outfile, 'w') as fh:
            fh.write('label, pvalue, tscore\n')
            for r in results:
                fh.write("{},{},{}\n".format(r[0], r[1], r[2]))

    def _tstats(self, wt_data, mut_data):

        results = []

        for label, wt_values in wt_data.iteritems():
            mut_values = mut_data.get(label)
            if not mut_values:
                results.append((label, 'no mutant data', 'no mutant data'))
                continue

            tscore, pvalue = stats.ttest_ind(np.array(
                mut_values, dtype=np.float),
                np.array(wt_values, dtype=np.float))

            results.append((label, pvalue, tscore))

        # Sort results. Lowest pval first
        sorted_results = reversed(sorted(results, key=lambda pval: pval[2]))
        return sorted_results




    def _get_label_vols(self, file_path):

        data = de
        with open(file_path, 'r') as csvfile:
            reader = csv.reader(csvfile)
            header = reader.next()
            for row in reader:
                for label_name, vol_calculation in zip(header[1:], row[1:], ):  # Skip the vol name column
                    data[label_name].append(vol_calculation)
        return data








