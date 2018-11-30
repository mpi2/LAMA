"""
Here we have some classes that hld all the information about a given statistical test
- input and output paths
- References to the input and putput data
"""

from pathlib import Path
from typing import Dict
from collections import defaultdict
import subprocess as sub
import numpy as np
import addict
import tempfile
from logzero import logger as logging

from lama import common
from lama.stats.standard_stats.data_loaders import InputData

RSCRIPT_FDR = common.lama_root_dir / 'stats' / 'rscripts' / 'r_padjust.R'

R_CHUNK_SIZE = 5000000


class Stats():
    specimen_results: addict.Dict

    def __init__(self,
                 input_: InputData,
                 stats_type: str
                 ):
        """

        Parameters
        ----------
        input_
            The data and associated metadata
        stats_type:
            jacobian, intensity etc
        """
        self.input_ = input_
        self.stats_type_ = stats_type
        self.stats_runner = None

        # The results will be store in these attributes
        self.line_qvals = None
        self.line_tstats = None
        self.specimen_results = None

    @staticmethod
    def factory(type_):
        if type_ == 'intensity':
            return Intensity
        elif type_ == 'jacobians':
            return Jacobians
        elif type_ == 'organ_volume':
            return OrganVolume

    def run_stats(self):
        """"""
        # These contain the chunked stats results
        line_level_pvals = []
        line_level_tvals = []

        # These contain the specimen-level statstics
        specimen_pvals = defaultdict(list)
        specimen_tstats = defaultdict(list)

        info = self.input_.info

        for data_chunk in self.input_.chunks(R_CHUNK_SIZE):

            current_chunk_size = data_chunk.shape[1]  # Not all chunks wil be same size

            p_all, t_all = self.stats_runner(data_chunk, info)

            # Convert all NANs in the pvalues to 1.0. Need to check that this is appropriate
            p_all[np.isnan(p_all)] = 1.0

            # Convert NANs to 0. We get NAN when for eg. all input values are 0
            t_all[np.isnan(t_all)] = 0.0

            # The first chunk of data will be from the line-level call
            p_line = p_all[:current_chunk_size]
            t_line = t_all[:current_chunk_size]

            line_level_pvals.append(p_line)
            line_level_tvals.append(t_line)

            # Get the specimen-level statistics
            mut_ids = self.input_.mutant_ids()

            for r, id_ in enumerate(mut_ids):
                start = current_chunk_size * (r + 1)
                end = current_chunk_size * (r + 2)
                t = t_all[start:end]
                p = p_all[start:end]
                specimen_tstats[id_].append(t)
                specimen_pvals[id_].append(p)

        # Stack the results chunks column-wise to get back to orginal shape
        line_pvals_array = np.hstack(line_level_pvals)
        line_tvals_array = np.hstack(line_level_tvals)

        self.line_qvals = fdr(line_pvals_array)

        self.line_tstats = line_tvals_array

        # Join up the results chunks for the specimen-level analysis. Do FDR correction on the pvalues
        self.specimen_results = addict.Dict()

        for id_, pvals in list(specimen_pvals.items()):
            p_ = np.array(pvals).ravel()
            q_ = fdr(p_)
            t_ = np.array(specimen_tstats[id_]).ravel()
            self.specimen_results[id_]['histogram'] = np.histogram(p_, bins=100)[0]
            self.specimen_results[id_]['q'] = q_
            self.specimen_results[id_]['t'] = t_
            self.specimen_results[id_]['p'] = p_



class Intensity(Stats):
    def __init__(self, *args):
        super().__init__(*args)


class Jacobians(Stats):
    def __init__(self, *args):
        super().__init__(*args)


class OrganVolume(Stats):
    def __init__(self, *args):
        super().__init__(*args)


def fdr(pvals):
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

    qval_outfile = tempfile.NamedTemporaryFile().name
    pval_file = tempfile.NamedTemporaryFile().name

    try:
        pvals.tofile(pval_file)
    except IOError:
        pass

    cmdFDR = ['Rscript',
              str(RSCRIPT_FDR),
              pval_file,
              str(pvals.shape[0]),
              qval_outfile
              ]

    try:
        sub.check_output(cmdFDR)
    except sub.CalledProcessError as e:
        logging.warn(f"R FDR calculation failed: {e}")
        raise

    return np.fromfile(qval_outfile, dtype=np.float64).astype(np.float32)

