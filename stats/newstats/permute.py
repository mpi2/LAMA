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
from pathlib import Path
import sys
sys.path.insert(0, Path(__file__).absolute() / '..')
from linear_model import lm_r


home = expanduser('~')



def get_mutant_n(dir_):
    """

    Parameters
    ----------
    dir_

    Yields
    -------
    int
        mutant n per line

    """
    mutant_ns = []
    for line_dir in [join(dir_, x) for x in os.listdir(dir_)]:
        if not os.path.isdir(line_dir): continue
        organ_vol_csv = join(line_dir, 'output', 'raw_label_sizes.csv')

        if not os.path.isfile(organ_vol_csv):
            print(f"can't find {organ_vol_csv}")
            continue

        df = pd.read_csv(organ_vol_csv, index_col=0)

        df = df[~df.index.str.contains('wt', case=False)]

        mutant_ns.append(len(df))
    return mutant_ns



def permuter(wt_csv:Path,
             wt_staging_csv:Path,
             mutant_dir: Path,
             out_path: Path,
             num_perm:int,
             plot_dir: Union[None, Path],
             specimen_level: bool=False,
             boxcox: bool=False,
             log_dependent: bool=False,
             log_staging: bool=False):
    """
    Main function loop over the mutant n and create syntheitc mutants and store the pvalues

    Parameters
    ----------
    wt_csv: csv containing the
    wt_staging_csv: csv containing the staging 
    mutant_dir
    out_path
    num_perm
    plot_dir
    specimen_level
    boxcox
    log_dependent
    log_staging

    Returns
    -------

    """
    # Get the crown-rump length for the baseliens
    df_crl = pd.read_csv(wt_staging_csv, index_col=0)

    df_crl.rename(columns={'value': 'crl'}, inplace=True)
    df_crl.index = df_crl.index.astype(str)

    if log_staging:
        print('logging staging metric')
        df_crl = np.log(df_crl)

    # Get the baseline data (organ volumes)
    df_wt = pd.read_csv(wt_csv, index_col=0)
    df_wt.index = df_wt.index.astype(str)

    df = df_wt.merge(right=df_crl, left_index=True, right_index=True)

    # We need to convert all the column names to not begin with a digit as LM does not like this
    df.columns = [f'x{x}'if x.isdigit() else x for x in df.columns]


    # Drop any columns that are all 0. These are the gaps in the label map caused by me merging labels
    df = df.loc[:, (df != 0).any(axis=0)]

    if log_dependent:
        print('logging dependent variable')
        df = np.log(df)

    label_header = df.drop(['crl'], axis='columns').columns

    p_results =[]

    # keep a list of sets of synthetic mutants, only run a set once
    synthetics_sets_done = []

    random.seed(999)

    if specimen_level:
        # Loop over the baselines, each time relablelling 1 baseline as synthetic mutant
        for n in range(len(df)):
            print(f'specimen-level {n}')
            df['genotype'] = 'wt'
            df.ix[[n], 'genotype'] = 'synth_hom'

            p = lm_r(df, plot_dir, boxcox)
            p_results.append(p)
    else:
        perms_done = 0
        for _ in range(num_perm):

            for n in get_mutant_n(mutant_dir): # muant lines
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

    pvalues_df = pd.DataFrame.from_records(p_results, columns=label_header)
    # Get rid of the x in the headers
    pvalues_df.columns = [x.strip('x') for x in pvalues_df.columns]
    pvalues_df.to_csv(out_path)




