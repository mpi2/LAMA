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


home = expanduser('~')




def plotter(qvalues_csv, out_dir):
    """
    This function is just here to be run after the permuting to make p-value distribution plots
    """

    q_df = pd.read_csv(qvalues_csv, index_col=0)

    for label in q_df:
        outpath = join(out_dir, f'{label}.png')
        ax = sns.distplot(q_df[label], bins=50)
        ax.set(xlim=(0, 1))
        plt.savefig(outpath)
        plt.close()



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


def get_all_files(root_dir, endswith):
    for root, dirs, files in os.walk(root_dir):
        for name in files:
            if name.endswith(endswith):
                yield os.path.join(root, name)


def prepare_mutant_data(mutant_dir, organ_info, out_path, stats_folder_name, specimen_level=False):
    """
    Permuter takes all the baseline from all specimens in csv files.
    This function just aggregates the mutant p value data from individual files so data so permuter can use it

    format of new dataframe should be specimen_id(index), label_num, label_num
    """
    organ_vol_pvalue_dfs = []

    def add_csv(path, first):
        if 'het' in path.lower():  # This removes lline
            print(path)
            return

        # extract p-values
        df = pd.read_csv(path, index_col=0)

        organ_vol_pvalue_dfs.append(df[['p']].T)  # .T to make organs as columns

        cols = df[['p']].T.columns
        if first:
            test = set(cols)
            first = False
        else:
            if test != set(cols):
                raise ValueError(f'{organ_vol_csv} has a different index than the rest(shoud be organ names)')

    first = True

    for line_dir in [join(mutant_dir, x) for x in os.listdir(mutant_dir)]:
        if not os.path.isdir(line_dir): continue

        if not specimen_level:  # Get pvalues for the line level results
            organ_vol_csv = join(line_dir, 'output', 'stats', stats_folder_name, 'inverted_organ_volumes_LinearModel_FDR5%.csv')
            if not os.path.isfile(organ_vol_csv):
                print(f"Can't find {organ_vol_csv}")
                continue
            add_csv(organ_vol_csv, first)
        else:
            speciemn_level_dir = join(line_dir, 'output', 'stats', stats_folder_name, 'specimen_calls')
            for spec_csv in get_all_files(speciemn_level_dir, '.csv'):
                add_csv(spec_csv, first)

    # merge into one dataframe
    df_result = pd.concat(organ_vol_pvalue_dfs, axis=0)

    # Read in label info csv
    labels_df = pd.read_csv(organ_info)

    # Replace label names with label numbers in df_result
    # cols = [labels_df[labels_df.label_name == x]['label'].values[0] for x in df_result.columns ]

    # cols = []
    # for label in df_result.columns:
    #     try:
    #         label_num = labels_df[labels_df.label_name == name]['label'].values[0]
    #     except IndexError as e:
    #         print(e)
    #     cols.append(label_num)



    # df_result.columns = cols

    df_result = df_result.sort_index(axis=1)


    df_result.to_csv(out_path)

    # The previous WT result used label number not name, so convert


def permuter(wt_csv, wt_crl_csv, mutant_dir, out_path, num_perm, plot_dir=None, specimen_level=False, boxcox=False, log_dependent=False, log_staging=False):
    """
    Main function loop over the mutant n and create syntheitc mutants and store the pvalues

    Parameters
    ----------
    mutant_dir: str
    wt_csv: str
        path to wt_organ vol csv

    Returns
    -------

    """
    # Get the crown-rump length for the baseliens
    df_crl = pd.read_csv(wt_crl_csv, index_col=0)

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




