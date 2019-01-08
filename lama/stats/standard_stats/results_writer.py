"""
This module does some post-processing of the stats results
 and writes out the results to file

 TODO: Remove some duplicated code
"""

from pathlib import Path
from typing import Tuple, List

import logzero
from logzero import logger as logging
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from lama.common import write_array, date_dhm
from lama.stats.standard_stats.stats_objects import Stats
from lama.stats.standard_stats.data_loaders import InputData

MINMAX_TSCORE = 50
FDR_CUTOFF = 0.05


class ResultsWriter:
    def __init__(self,
                 results: Stats,
                 mask: np.ndarray,
                 root_out_dir: Path,
                 stats_name: str,
                 label_info_path: Path):
        """
        TODO: map organ names back onto results
        Parameters
        ----------
        results
            The object containing all the stats results
        mask
            Mask. Not needed for organ volumes
        root_out_dir
            The root directory to create a subdirectory in to store output.
        stats_name
            The name of the type of analysis (eg intensity)
        label_info_file
            Label map information

        Returns
        -------

        """
        self.label_info_path = label_info_path
        self.results = results
        self.mask = mask
        self.shape = results.input_.shape
        self.stats_name = stats_name
        self.line = results.input_.line

        # Write out the line-level results
        line_tstats = results.line_tstats
        line_qvals = results.line_qvals

        self.out_dir = root_out_dir / self.line/ stats_name
        self.out_dir.mkdir(parents=True, exist_ok=True)

        self._write(line_tstats, line_qvals, self.line)

        # For out specimen-level results
        for spec_id, spec_res in results.specimen_results.items():

            spec_t = spec_res['t']
            spec_q = spec_res['q']
            self._write(spec_t, spec_q, spec_id)

        # self.log(self.out_dir, 'Organ_volume stats', results.input_)

    @staticmethod
    def factory(data_type):

        return {'jacobians': VoxelWriter,
                'intensity': VoxelWriter,
                'organ_volumes': OrganVolumeWriter
                }[data_type]

    def log(self):
        logfile = self.out_diroutdir / f'{date_dhm()}_stats.log'
        logzero.logfile(logfile)

        logging.info(f'Doing {self.stats_name}')
        wt_paths = '\n'.join([str(x) for x in self.results.input_.paths[0]])
        mut_patsh = '\n'.join([str(x) for x in self.results.input_.paths[1]])
        logging.info(f'wild type inputs: {wt_paths}')
        logging.info(f'Mutant inputs: {mut_patsh}')

    def _write(self, t_stats, qvals, name):
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
        super().__init__(*args)

    def _write(self, t_stats, qvals, name):

        filtered_tstats = result_cutoff_filter(t_stats, qvals)
        line_result = self.rebuild_array(filtered_tstats, self.shape, self.mask)

        out_path = self.out_dir/ f'{name}_{self.stats_name}.nrrd'
        write_array(line_result, out_path)

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

    def _write(self, t_stats, qvals, name):
        # write_csv(self.line_tstats, self.line_qvals, line_out_path, list(results.input_.data.columns), label_info)
        out_path = self.out_dir / f'{name}_{self.stats_name}.csv'
        df = pd.DataFrame.from_dict(dict(t=t_stats, q=qvals))

        label_info = pd.read_csv(self.label_info_path)

        # Merge the results from each label to the label info
        # The labels are stored in the InputData
        labels = list(self.results.input_.data.columns)
        df.index = labels
        df.index = df.index.astype(np.int64)
        df = df.merge(right=label_info, right_on='label', left_index=True)

        df['significant_bh_q_5%'] = df['q'] < 0.05
        df.to_csv(out_path)


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
        t[mask] = 0

    return t


def pvalue_fdr_plot(pvals, outdir: Path):
    """
    Write out a fdr correction plot.

    Got the idea from: https://www.unc.edu/courses/2007spring/biol/145/001/docs/lectures/Nov12.html
    """
    # Make pvalue plot

    line_fdr_fig = outdir / 'fdr_correction.png'
    sorted_p = sorted(list(pvals))
    x = np.array(list(range(len(sorted_p)))) / float(len(sorted_p))  # k_m = rank/num pvalues
    sns.scatterplot(x=x, y=sorted_p, sizes=(1,))
    plt.plot([0, 1], [0, 0.05])
    plt.xlabel('k/m')
    plt.ylabel('p-value')
    plt.savefig(line_fdr_fig)
    plt.close()
