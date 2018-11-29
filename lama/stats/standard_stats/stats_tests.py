"""
This module does the statsistical analysis
"""

import subprocess as sub
import tempfile
import numpy as np
from collections import defaultdict

import pandas as pd
from logzero import logger as logging

from lama.stats.standard_stats.data_loaders import InputData
from lama.stats.standard_stats.linear_model import lm_r

R_CHUNK_SIZE = 5000000

class OutputData:
    def __init__(self):
        pass



def lm_r(input_: InputData):
    """
    Given an InputData instance regress the value at each columns against genotype and optionally staging metric
    """
    # Make temp files for storing the output of the line levl LM results
    line_level_pval_out_file = tempfile.NamedTemporaryFile().name
    line_level_tstat_out_file = tempfile.NamedTemporaryFile().name

    # logging.info('using formula {}'.format(self.formula))
    # print('num chunks', num_chunks)

    # Loop over the data in chunks
    # np.array_split provides a split view on the array so does not increase memory
    # The result will be a bunch of arrays split across the second dimension
    # chunked_data = np.array_split(data, num_chunks, axis=1)


    # These contain the chunked stats results
    line_level_pvals = []
    line_level_tvals = []

    # These contain the specimen-level statstics
    specimen_pvals = defaultdict(list)
    specimen_tstats = defaultdict(list)

    # Load in the groups file so we can get the specimen names
    # groups_df = pd.read_csv(self.groups)
    # mutants_df = groups_df[groups_df.genotype == 'mutant']

    i = 0  # enumerate!
    voxel_file = tempfile.NamedTemporaryFile().name
    # voxel_file = '/home/neil/Desktop/t/test_stats/test_voxel_file_67_68'

    # plot_out_dir = Path(self.outdir) / 'lm_plots'
    # plot_out_dir.mkdir(exist_ok=True)

    groups = input_.info
    for data_chunk in input_.chunks(R_CHUNK_SIZE):

        current_chunk_size = data_chunk.shape[1]  # Not all chunks wil be same size
        print('chunk: {}'.format(i))


        lm_r(data_chunk, groups)

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



def fdr_bh(pvals):
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
