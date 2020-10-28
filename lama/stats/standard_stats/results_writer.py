"""
This module does some post-processing of the stats results
 and writes out the results to file

ResultsWriter uis subclassed for each data type.
Currently just the *_write* method is overridden in subclasses
 which take into account the differences in output between voxel based data (3D volume) and organ volume results (CSV file)


"""

from pathlib import Path
from typing import Tuple

import logzero
from logzero import logger as logging
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from lama.common import write_array
from lama.stats.standard_stats.stats_objects import Stats

MINMAX_TSCORE = 50
FDR_CUTOFF = 0.05

# 041219
# The stats files lose the header information and are written with an incorrect lps header wothout flippin gthe spaces
# For now override this behaviour by adding RAS header. But this assumes the input is in RAS. Fix ASAP


class ResultsWriter:
    def __init__(self,
                 results: Stats,
                 mask: np.ndarray,
                 out_dir: Path,
                 stats_name: str,
                 label_map: np.ndarray,
                 label_info_path: Path):
        """
        TODO: map organ names back onto results
        Parameters
        ----------
        results
            The object containing all the stats results
        mask
            Mask. Not needed for organ volumes
        out_dir
            The root directory to create a subdirectory in to store output.
        stats_name
            The name of the type of analysis (eg intensity)
        label_map
            for creating filtered labelmap overlays
        label_info_path
            Label map information

        Returns
        -------

        """
        self.label_info_path = label_info_path
        self.label_map = label_map
        self.out_dir = out_dir
        self.results = results
        self.mask = mask
        self.shape = results.input_.shape
        self.stats_name = stats_name
        self.line = results.input_.line

        # Write out the line-level results
        line_tstats = results.line_tstats
        line_qvals = results.line_qvals
        line_pvals = results.line_pvalues # Need to get thse into organ volumes

        logging.info('Writing line-level results...')
        
        line_threshold_file = self.out_dir / f'Qvals_{stats_name}_{self.line}.csv'
        write_threshold_file(line_qvals, line_tstats, line_threshold_file)

        # this is for the lama_stats to know where the heatmaps for inversion are
        self.line_heatmap = self._write(line_tstats, line_pvals, line_qvals, self.out_dir, self.line)  # Bodge. Change!
        # 140620: I think this is where it's dying
        # pvalue_fdr_plot(results.line_pvalues, results.line_qvals, out_dir)

        logging.info('Finished writing line-level results')

        logging.info('Writing specimen-level results...')
        
        specimen_out_dir = out_dir / 'specimen-level'
        specimen_out_dir.mkdir(exist_ok=True)

        # For specimen-level results
        for spec_id, spec_res in results.specimen_results.items():
            spec_threshold_file = specimen_out_dir / f'Qvals_{stats_name}_{spec_id}.csv'
            spec_t = spec_res['t']
            spec_q = spec_res['q']
            spec_p = spec_res['p']
            write_threshold_file(spec_q, spec_t, spec_threshold_file)
            self._write(spec_t, spec_p, spec_q, specimen_out_dir, spec_id)

        # self.log(self.out_dir, 'Organ_volume stats', results.input_)
        logging.info('Finished writing specimen-level results')


    @staticmethod
    def factory(data_type):

        return {'jacobians': VoxelWriter,
                'intensity': VoxelWriter,
                'organ_volumes': OrganVolumeWriter
                }[data_type]

    def _write(self):
        """
        Write the results to file
        """
        raise NotImplementedError


class VoxelWriter(ResultsWriter):
    def __init__(self, *args):
        """
         Write the line and specimen-level results.

         Remove any Nans
         Threshold the t-statstistics based on q-value
         Write nrrds to file.

         Parameters
         ----------
         results
         out_dir
             The root directory to put the results in
         stats_name
             An the stats type
         label_info:
             Not currently used
         """
        self.line_heatmap = None
        super().__init__(*args)

    def _write(self, t_stats, pvals, qvals, outdir, name):
        filtered_tstats = result_cutoff_filter(t_stats, qvals)
        filtered_result = self.rebuild_array(filtered_tstats, self.shape, self.mask)
        unfiltered_result = self.rebuild_array(t_stats, self.shape, self.mask)

        heatmap_path = outdir / f'{name}_{self.stats_name}_t_fdr5.nrrd'
        heatmap_path_unfiltered = outdir / f'{name}_{self.stats_name}_t.nrrd'

        # Write qval-filtered t-stats
        write_array(filtered_result, heatmap_path, ras=True)

        # Write raw t-stats
        write_array(unfiltered_result, heatmap_path_unfiltered, ras=True)

        return heatmap_path


    @staticmethod
    def rebuild_array(array: np.ndarray, shape: Tuple, mask: np.ndarray) -> np.ndarray:
        """
        The stats pipeline uses masked data throughout to save on resources
        This function rebuilds the output files to the orginal 3D sizes of the input volumes

        Parameters
        ----------
        array
            1d masked array to rebuild
        shape
            shape of input volume
        mask
            3D mask

        Returns
        -------
        3d rebuilt array

        """

        array[array > MINMAX_TSCORE] = MINMAX_TSCORE
        array[array < -MINMAX_TSCORE] = - MINMAX_TSCORE

        full_output = np.zeros(shape)
        full_output[mask != False] = array
        return full_output.reshape(shape)


class OrganVolumeWriter(ResultsWriter):
    def __init__(self, *args):
        super().__init__(*args)
        self.line_heatmap = None

        # Expose the results for clustering
        self.organ_volume_results: pd.DataFrame = None#????

    def _write(self, t_stats, pvals, qvals, out_dir, name):
        # write_csv(self.line_tstats, self.line_qvals, line_out_path, list(results.input_.data.columns), label_info)
        out_path = out_dir / f'{name}_{self.stats_name}.csv'
        df = pd.DataFrame.from_dict(dict(t=t_stats, p=pvals, q=qvals))

        label_info = pd.read_csv(self.label_info_path)

        # Merge the results from each label to the label info
        # The labels are stored in the InputData
        labels = list(self.results.input_.data.columns)
        df.index = labels
        df.index = df.index.astype(np.int64)
        df = df.merge(right=label_info, right_on='label', left_index=True)

        df['significant_bh_q_5'] = df['q'] < 0.05
        df.sort_values(by='q', inplace=True)
        df.to_csv(out_path)

        hit_labels = df[df['significant_bh_q_5'] == True]['label']

        thresh_labels_out = out_dir / f'{name}_hit_organs.nrrd'
        # self._write_thresholded_label_map(self.label_map, hit_labels, thresh_labels_out)

    def _write_thresholded_label_map(self, label_map: np.ndarray, hits, out: Path):
        """
        Write a label map with only the 'hit' organs in it
        """
        if len(hits) > 0:
            # Make a copy as it may be being used elsewhere
            l = np.copy(label_map)
            # Clear any non-hits
            l[~np.isin(l, hits)] = 0

            write_array(l, out, ras=True)


def result_cutoff_filter(t: np.ndarray, q: np.ndarray) -> np.ndarray:
    """
    Convert to numpy arrays and set to zero any tscore that has a corresponding pvalue > 0.05

    Parameters
    ----------

    """
    if len(t) != len(q):
        raise ValueError
    else:
        mask = q > FDR_CUTOFF
        masked = np.copy(t)
        masked[mask] = 0

        return masked


def pvalue_fdr_plot(pvals, qvals, outdir: Path):
    """
    Write out a fdr correction plot.

    Got the idea from: https://www.unc.edu/courses/2007spring/biol/145/001/docs/lectures/Nov12.html
    """
    # Make p-value plot

    line_fdr_fig = outdir / 'fdr_correction.png'

    # debug.


    df = pd.DataFrame.from_dict({'p': pvals, 'q': qvals})
    df.sort_values('p', ascending=True)
    df['k/m'] = np.array(list(range(len(df)))) / float(len(df))  # k_m = rank/num pvalues

    df['significant'] = df.q <= 0.05

    ax = sns.scatterplot(x='k/m', y='p', hue='significant', data=df)

    sns.lineplot([0, 1], [0, 0.05])
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles, ["", "Significant", "Not significant"])
    plt.savefig(line_fdr_fig)
    plt.close()


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
