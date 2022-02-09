"""
The code to run linear models in R.

The current interface to R is to write binary files that R can read (numpy_to_dat). The reason r2py wasn't used is that
it used to be a pain to install. I imagine it's better now and using docker should improve things so adding
rp2y interface is on the todo list
"""

import subprocess as sub
import os
import struct
from pathlib import Path
import tempfile
from typing import Tuple
import shutil

from itertools import chain
from logzero import logger as logging

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from lama import common

LM_SCRIPT = str(common.lama_root_dir / 'stats' / 'rscripts' / 'lmFast.R')

# If debugging, don't delete the temp files used for communication with R so they can be used for R debugging.
DEBUGGING = True


def lm_r(data: np.ndarray, info: pd.DataFrame, plot_dir: Path = None, boxcox: bool = False, use_staging: bool = True,
         two_way: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit multiple linear models and get the resulting p-values

    Parameters
    ----------
    data
        columns: data points
        rows: specimens
    info
        columns:
            genotype, staging, line
        rows:
            specimens
    plot_dir
        where to optionally output lm plots (qq etc)
    boxcox
        whether to apply boxcox transformation to the dependent variable
    use_staging
        if true, uae staging as a fixed effect in the linear model

    Returns:
    -------
    pvalues for each label or voxel
    t-statistics for each label or voxel

    """

    input_binary_file = tempfile.NamedTemporaryFile().name
    line_level_pval_out_file = tempfile.NamedTemporaryFile().name
    line_level_tstat_out_file = tempfile.NamedTemporaryFile().name
    groups_file = tempfile.NamedTemporaryFile().name

    # create groups file
    if use_staging and two_way:
        groups = info[['genotype', 'treatment', 'staging']]
        formula = 'genotype,treatment,staging'

    elif use_staging:
        groups = info[['genotype', 'staging']]
        formula = 'genotype,staging'

    else:
        groups = info[['genotype']]
        formula = 'genotype'

    groups.index.name = 'volume_id'
    groups.to_csv(groups_file)

    _numpy_to_dat(data, input_binary_file)

    cmd = ['Rscript',
           LM_SCRIPT,
           input_binary_file,
           groups_file,
           line_level_pval_out_file,
           line_level_tstat_out_file,
           formula,
           str(boxcox).upper(),  # bool to string for R
           ''  # No plots needed for permutation testing
           ]

    try:
        sub.check_output(cmd)
        logging.info('R linear model suceeded')
    except sub.CalledProcessError as e:
        msg = "R linear model failed: {}".format(e)
        logging.exception(msg)
        raise RuntimeError(msg)

    # Read in the pvalue and t-statistic results.
    # The start of the binary file will contain values from the line level call
    # the specimen-level calls are appended onto this and need to be split accordingly.
    try:

        if two_way:
            # if you're doing a two-way, you need to get three arrays (i.e. genotype, treatment and interaction)
            # converted to np.ndarray with three dimensions 

            p_all = np.array(
                [np.fromfile((line_level_pval_out_file + '_genotype'), dtype=np.float64).astype(np.float32),
                 np.fromfile((line_level_pval_out_file + '_treatment'), dtype=np.float64).astype(np.float32),
                 np.fromfile((line_level_pval_out_file + '_interaction'), dtype=np.float64).astype(np.float32)])

            t_all = np.array(
                [np.fromfile((line_level_tstat_out_file + '_genotype'), dtype=np.float64).astype(np.float32),
                 np.fromfile((line_level_tstat_out_file + '_treatment'), dtype=np.float64).astype(np.float32),
                 np.fromfile((line_level_tstat_out_file + '_interaction'), dtype=np.float64).astype(np.float32)])
        else:
            p_all = np.fromfile(line_level_pval_out_file, dtype=np.float64).astype(np.float32)
            t_all = np.fromfile(line_level_tstat_out_file, dtype=np.float64).astype(np.float32)

    except FileNotFoundError as e:
        print(f'Linear model file from R not found {e}')
        raise FileNotFoundError('Cannot find LM output'.format(e))

    if not DEBUGGING:
        os.remove(input_binary_file)
        os.remove(groups_file)

        if two_way:
            os.remove(line_level_pval_out_file + '_genotype')
            os.remove(line_level_pval_out_file + '_treatment')
            os.remove(line_level_pval_out_file + '_interaction')
            os.remove(line_level_tstat_out_file + '_genotype')
            os.remove(line_level_tstat_out_file + '_treatment')
            os.remove(line_level_tstat_out_file + '_interaction')
        else:
            os.remove(line_level_pval_out_file)
            os.remove(line_level_tstat_out_file)
    return p_all, t_all


def _numpy_to_dat(mat: np.ndarray, outfile: str):
    """
    Convert a numpy array to a binary file for reading in by R

    Parameters
    ----------
    mat
        the data to be send to r
    outfile
        the tem file name to store the binary file


    """
    # mat = mat.as_matrix()
    # create a binary file
    with open(outfile, 'wb') as binfile:
        # and write out two integers with the row and column dimension

        header = struct.pack('2I', mat.shape[0], mat.shape[1])
        binfile.write(header)
        # then loop over columns and write each
        for i in range(mat.shape[1]):
            data = struct.pack('%id' % mat.shape[0], *mat[:, i])
            binfile.write(data)


def lm_sm(data: np.ndarray, info: pd.DataFrame, plot_dir: Path = None,
          boxcox: bool = False, use_staging: bool = True, two_way: bool = False):
    """

    Parameters
    ----------
    two_way
    data
    info
    plot_dir
    boxcox
    use_staging

    Notes
    -----
    If a label column is set to all 0, it means a line has all the mutants qc's and it's not for analysis.

    Returns
    -------

    """

    pvals = []
    tvals = []

    # We need to add some non-digit before column names or patsy has a fit
    d = pd.DataFrame(data, index=info.index, columns=[f'x{x}' for x in range(data.shape[1])])
    df = pd.concat([d, info], axis=1)  # Data will be given numberic labels

    for col in range(data.shape[1]):

        if not df[f'x{col}'].any():
            p = np.nan
            t = np.nan
        else:
            if two_way:
                # two way model - if its just the geno or treat comparison; the one-factor col will
                # be ignored
                # for some simulations smf is being a cunt and is returning the text in pvalues.

                fit = smf.ols(formula=f'x{col} ~ genotype * treatment + staging', data=df, missing='drop').fit()

                # get all pvals except intercept and staging

                p = fit.pvalues[~fit.pvalues.index.isin(['Intercept', 'staging'])]
                t = fit.tvalues[~fit.tvalues.index.isin(['Intercept', 'staging'])]
            else:
                fit = smf.ols(formula=f'x{col} ~ genotype + staging', data=df, missing='drop').fit()
                p = fit.pvalues['genotype[T.wt]']
                t = fit.tvalues['genotype[T.wt]']
        pvals.append(p)
        tvals.append(t)

    # print(dir(fit))

    p_all = np.array(pvals)
    # debugging which may not matter

    # weird output from smf.ols.fit.pvalues
    if len(np.shape(p_all)) == 1:
        # Coerces all the data to be the right shape

        p_all = np.reshape(p_all, (187, 1))

    # print(type(p_all[0][0]), p_all[0][0])
    if isinstance(p_all[0][0], pd.Series):
        # coerces the series in p_all to np.float
        p_all = [np.float64(my_p) for my_pvals in p_all for my_p in my_pvals]

    # now reshape properly
    if len(np.shape(p_all)) == 1:
        # Coerces all the data to be the right shape (i.e 187x1 or 187x3)
        if isinstance(p_all[0], np.ndarray):
            p_all = [np.array(3*[np.nan]) if np.isnan(org).any() else org for org in p_all]
            p_all = np.vstack(p_all)
        else:
            p_all = np.reshape(p_all, (187, 1))

    t_all = np.negative(np.array(tvals))  # The tvaue for genotype[T.mut] is what we want

    return p_all, t_all
