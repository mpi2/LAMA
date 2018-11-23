#!/usr/bin/env python3

"""
06/08/18

* Determine a null distributiom in order to set an appropriate p-value threshold
* Try first on organ volume results

Steps
=====

For each mutant line
    get n number of homs
    relabel n wildtypes as mutants
    run organ volume LM   organ_volume ~ genotype + 'crown-rump length'


See labbook with this date for futther details

"""
from typing import Dict
import os
from os.path import expanduser, join
from typing import Union
from collections import defaultdict
import statsmodels.formula.api as smf
import pandas as pd
import random
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess as sub
import struct
from logzero import logger as logging
from pathlib import Path
import sys
sys.path.insert(0, Path(__file__).absolute() / '..')
from linear_model import lm_r


home = expanduser('~')


def permuter(wt_organ_vol_csv: Path,
             wt_staging_csv:Path,
             mut_staging_csv: Path,
             out_path: Path,
             num_perm:int,
             plot_dir: Union[None, Path]=None,
             specimen_level: bool=False,
             boxcox: bool=False,
             log_dependent: bool=False,
             log_staging: bool=False):
    """
    Main function loop over the mutant n and create syntheitc mutants and store the pvalues

    Parameters
    ----------
    wt_organ_vol_csv
        csv containing the organ volumes
        columns: vol, label numbers... line
    wt_staging_csv
        csv containing the staging
        columns: vol, value(the staging metric), line
    mut_staging_csv
        csv listing all mutant specimens
        Use the 'line' column to get the mutant n number for each line

    out_path
    num_perm
    plot_dir
    specimen_level
    boxcox
    log_dependent
    log_staging

    Returns
    -------

    Notes
    -----
    TODO: Remove crl and using 'staging instead'

    """
    # Get the crown-rump length for the baseliens
    df_crl = pd.read_csv(wt_staging_csv, index_col=0)

    df_crl.rename(columns={'value': 'crl'}, inplace=True)
    df_crl.index = df_crl.index.astype(str)

    if log_staging:
        print('logging staging metric')
        df_crl = np.log(df_crl)

    # Get the baseline data (organ volumes)
    df_wt = pd.read_csv(wt_organ_vol_csv, index_col=0)

    # remove the line columns as it's not needed as it will just be 'baseline'
    df_wt = df_wt.drop(columns=['line'])
    df_crl = df_crl.drop(columns=['line'])

    # Drop any columns that are all 0. These are the gaps in the label map caused by me merging labels
    df_wt = df_wt.loc[:, (df_wt != 0).any(axis=0)]

    # For the statsmodels linear mode to work, column names cannot start with a digid. Prefix with 'x'
    df_wt.columns = [f'x{x}' if x.isdigit() else x for x in df_wt.columns]

    if log_dependent:
        logging.info('logging dependent variable')
        df_wt = np.log(df_wt)

    # Index needs to be a str for the merge to work
    df_wt.index = df_wt.index.astype(str)

    # Create a dataframe with organ volume columns + the staging column
    df = df_wt.merge(right=df_crl, left_index=True, right_index=True)

    df_mut = pd.read_csv(mut_staging_csv, index_col=0)

    # Get the line specimen n numbers. Keep the fits column
    line_specimen_counts = df_mut.groupby('line').count()
    line_specimen_counts = list(line_specimen_counts.iloc[:, 0])

    label_names = df.drop(['crl'], axis='columns').columns

    # Store p-value results. One tuple(len=num labels) per iteration
    p_results = []

    # keep a list of sets of synthetic mutants, only run a set once
    synthetics_sets_done = []

    random.seed(999)

    if specimen_level:
        # Create synthetic specimens by iteratively relabelling each basline as synthetic mutant
        for n in range(len(df)):
            print(f'specimen-level {n}')
            df['genotype'] = 'wt'
            df.ix[[n], 'genotype'] = 'synth_hom'

            p = lm_r(df, plot_dir, boxcox)
            p_results.append(p)
    else:
        # Create synthetic lines by iteratively relabelling n baselines as synthetic mutants
        # n is determined by sampling the number of homs in each mutant line
        perms_done = 0
        for _ in range(num_perm):

            for n in line_specimen_counts: # muant lines
                perms_done += 1
                print(perms_done)
                # Set all to wt genotype
                df['genotype'] = 'wt'

                # label n number of baselines as mutants
                synthetics_mut_indices = random.sample(range(0, len(df)), n)

                if synthetics_mut_indices in synthetics_sets_done:
                    continue

                synthetics_sets_done.append(synthetics_mut_indices)

                df.ix[synthetics_mut_indices, 'genotype'] = 'synth_hom'  # Why does iloc not work here?

                p = lm_r(df, plot_dir, boxcox) # returns p_values for all organs, 1 iteration
                p_results.append(p)

    pvalues_df = pd.DataFrame.from_records(p_results, columns=label_names)
    # Get rid of the x in the headers
    pvalues_df.columns = [x.strip('x') for x in pvalues_df.columns]
    pvalues_df.to_csv(out_path)




