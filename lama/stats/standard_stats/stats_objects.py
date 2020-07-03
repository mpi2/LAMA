"""
Here we have some classes that hold all the information about a given statistical test
- input and output paths
- References to the input and putput data

TODO
I thought that each data type might need it's own subclass, however the different DataLoaders and ResultWriters are enough
"""

from collections import defaultdict
import subprocess as sub
import tempfile
import os
from pathlib import Path

import numpy as np
import addict
from logzero import logger as logging
import pandas as pd


from lama import common
from lama.stats.standard_stats.data_loaders import LineData

RSCRIPT_FDR = common.lama_root_dir / 'stats' / 'rscripts' / 'r_padjust.R'


class Stats:
    specimen_results: addict.Dict

    def __init__(self,
                 input_: LineData,
                 stats_type: str,
                 use_staging: bool = True
                 ):
        """

        Parameters
        ----------
        input_
            The data and associated metadata
        stats_type:
            jacobian, intensity etc
        use_staging:
            add staging in linear model
        """
        self.input_ = input_
        self.stats_type_ = stats_type
        self.stats_runner = None
        self.use_staging = use_staging

        # The final results will be stored in these attributes
        self.line_qvals = None
        self.line_pvalues = None
        self.line_tstats = None
        self.specimen_results = None

    @staticmethod
    def factory(type_):
        if type_ == 'intensity':
            return Intensity
        elif type_ == 'jacobians':
            return Jacobians
        elif type_ == 'organ_volumes':
            return OrganVolume

    def run_stats(self):

        if self.use_staging:
            logging.info('Using genotype and staging in linear model')
        else:
            logging.info('Using only genotype in linear model')

        # These contain the chunked stats results
        line_level_pvals = []
        line_level_tvals = []

        # These will contain the specimen-level statstics
        specimen_pvals = defaultdict(list)
        specimen_tstats = defaultdict(list)

        info = self.input_.info

        num_chunks = self.input_.get_num_chunks(log=True)

        for i, data_chunk in enumerate(self.input_.chunks()):
            # Chunk the data and send sequentially to R to not use all the memory

            logging.info(f'Chunk {i + 1}/{num_chunks}')

            current_chunk_size = data_chunk.shape[1]  # Final chunk may not be same size

            p_all, t_all = self.stats_runner(data_chunk, info, use_staging=self.use_staging)

            # Convert all NANs in the pvalues to 1.0. Need to check that this is appropriate
            p_all[np.isnan(p_all)] = 1.0

            # Convert NANs to 0. We get NAN when for eg. all input values are 0
            t_all[np.isnan(t_all)] = 0.0

            # Each chunk of results has the line -level results at the start
            p_line = p_all[:current_chunk_size]
            t_line = t_all[:current_chunk_size]

            line_level_pvals.append(p_line)
            line_level_tvals.append(t_line)

            # Get the specimen-level statistics
            mut_ids = self.input_.mutant_ids()

            for spec_num, id_ in enumerate(mut_ids):
                # After the line level result, the specimen-level results are appended to the result chunk
                start = current_chunk_size * (spec_num + 1)
                end = current_chunk_size * (spec_num + 2)

                specimen_tstats[id_].append(t_all[start:end])
                specimen_pvals[id_].append(p_all[start:end])

        # Stack the results chunks column-wise to get back to orginal shape
        line_pvals_array = np.hstack(line_level_pvals)
        line_tvals_array = np.hstack(line_level_tvals)

        self.line_pvalues = line_pvals_array

        self.line_qvals = fdr(line_pvals_array)

        self.line_tstats = line_tvals_array

        # Join up the results chunks for the specimen-level analysis. Do FDR correction on the pvalues
        self.specimen_results = addict.Dict()

        try:
            for id_, p in list(specimen_pvals.items()):
                p = np.hstack(p)
                q = fdr(p)
                t = np.hstack(specimen_tstats[id_])
                self.specimen_results[id_]['histogram'] = np.histogram(p, bins=100)[0]
                self.specimen_results[id_]['q'] = q
                self.specimen_results[id_]['t'] = t
                self.specimen_results[id_]['p'] = p
        except Exception as e:
            logging.info(p)


class Intensity(Stats):
    def __init__(self, *args):
        super().__init__(*args)


class Jacobians(Stats):
    def __init__(self, *args):
        super().__init__(*args)


class OrganVolume(Stats):
    def __init__(self, *args):
        super().__init__(*args)


def fdr(pvals: np.ndarray) -> np.ndarray:
    """
    Use R for FDR correction

    ----------
    pvals: The p-values to be corrected

    Returns
    -------
    The corrected q-values
    """

    qval_outfile = tempfile.NamedTemporaryFile().name
    pval_file = tempfile.NamedTemporaryFile().name

    try:
        pvals.tofile(pval_file)
    except IOError:
        logging.error(f'Cannot save p-value file {pval_file}')


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

    result = np.fromfile(qval_outfile, dtype=np.float64).astype(np.float32)

    os.remove(qval_outfile)
    os.remove(pval_file)

    return result



