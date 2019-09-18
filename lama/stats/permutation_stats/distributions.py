#!/usr/bin/env python3

"""
06/08/18

* Determine a null distributiom in order to set an appropriate p-value threshold


Steps
=====

For each mutant line
    get n number of homs
    relabel n wildtypes as mutants
    run organ volume LM   organ_volume ~ genotype + 'staging metric'



"""

from os.path import expanduser
from typing import Union, Tuple, List
import random
from pathlib import Path

import pandas as pd
from scipy.special import comb
from logzero import logger as logging

from lama.stats.standard_stats.linear_model import lm_r

home = expanduser('~')


def null(input_data: pd.DataFrame,
         num_perm: int,
         plot_dir: Union[None, Path] = None,
         boxcox: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Generate null distributions for line and specimen-level data

    Parameters
    ----------
    input_data
        columns
            staging
            line
            then multiple columns each one a label (organ)
        Baselines must be labelled 'baseline' in the line column

    num_perm
        number of permutations
    plot_dir
        Where to store the optional lm plots
    boxcox
        If true, box cox transform dependent variable

    Returns
    -------
    line-level null distribution
    specimen-level null distribution

    Notes
    -----
    Labels must not start with a digit as R will throw a wobbly
    """
    random.seed(999)

    # Use the generic staging label from now on
    input_data.rename(columns={'crl': 'staging', 'volume': 'staging'}, inplace=True)

    label_names = input_data.drop(['staging', 'line'], axis='columns').columns

    # Store p-value results. One tuple(len=num labels) per iteration
    line_p = []
    line_t = []
    spec_p = []
    spec_t = []

    # Create synthetic specimens by iteratively relabelling each baseline as synthetic mutant
    baselines = input_data[input_data['line'] == 'baseline']

    # Get the line specimen n numbers. Keep the first column
    line_specimen_counts = input_data[input_data['line'] != 'baseline'].groupby('line').count()
    line_specimen_counts = list(line_specimen_counts.iloc[:, 0])

    # Split data into a numpy array of raw data and dataframe for staging and genotype fopr the LM code
    data = baselines.drop(columns=['staging', 'line']).values
    info = baselines[['staging', 'line']]

    # Get the specimen-level null distribution
    # TODO: Why does the LM not try to return a specimen-level p value as well?
    for index, _ in info.iterrows():
        info['genotype'] = 'wt'                     # Set all genotypes to WT
        info.ix[[index], 'genotype'] = 'synth_hom'  # Set the ith baseline to synth hom

        # Get a p-value for each organ
        p, t = lm_r(data, info)

        # Check that there are equal amounts of p-value sthan there are data points
        if len(p) != data.shape[1]:
            raise ValueError(f'The length of p-values results: {data.shape[1]} does not match the length of the input data: {len(p)}')

        spec_p.append(p)
        spec_t.append(t)

    # Line-level null distribution
    # Create synthetic lines by iteratively relabelling n baselines as synthetic mutants
    # n is determined by sampling the number of homs in each mutant line

    # keep a list of sets of synthetic mutants, only run a set once
    synthetics_sets_done = []

    perms_done = 1

    while perms_done < num_perm:

        for n in line_specimen_counts:  # mutant lines

            perms_done += 1

            if perms_done == num_perm:
                break

            print(f'permutation: {perms_done}')

            if not _label_synthetic_mutants(info, n, synthetics_sets_done):
                continue

            p, t = lm_r(data, info)  # returns p_values for all organs, 1 iteration

            if len(p) != data.shape[1]:
                raise ValueError(
                    f'The length of p-values results: {data.shape[1]} does not match the length of the input data: {len(p)}')

            line_p.append(p)
            line_t.append(t)

    line_df = pd.DataFrame.from_records(line_p, columns=label_names)
    spec_df = pd.DataFrame.from_records(spec_p, columns=label_names)

    # Get rid of the x in the headers that were needed for R
    strip_x([line_df, spec_df])

    return line_df, spec_df


def _label_synthetic_mutants(info: pd.DataFrame, n: int, sets_done: List):
    """
    Given a dataframe of wild type data relabel n baselines as synthetic mutant in place.
    Keep track of combinations done in sets_done and do not duplicate

    Parameters
    ----------
    info
        dataframe with 'genotype' column
    n
        how many specimens to relabel
    sets_done
        Contains Sets of previously selected specimen IDs
    """

    # Set all to wt genotype
    info['genotype'] = 'wt'

    # label n number of baselines as mutants

    max_comb = int(comb(len(info), n))

    for i in range(max_comb):
        synthetics_mut_indices = random.sample(range(0, len(info)), n)
        i += 1
        if not set(synthetics_mut_indices) in sets_done:
            break

    if i > max_comb - 1:
        msg = f"""Cannot find unique combinations of wild type baselines to relabel as synthetic mutants
        With a baseline n  of {len(info)}\n. Choosing {n} synthetics. 
        Try increasing the number of baselines or reducing the number of permutations"""

        return False

    sets_done.append(set(synthetics_mut_indices))

    info.ix[synthetics_mut_indices, 'genotype'] = 'synth_hom'  # Why does iloc not work here?
    return True


def strip_x(dfs):
    for df in dfs:
        df.columns = [x.strip('x') for x in df.columns]


def alternative(input_data: pd.DataFrame,
                plot_dir: Union[None, Path] = None,
                boxcox: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Generate alterntive (mutant) distributions for line and pecimen-level data

    Parameters
    ----------
    input_data
    plot_dir
    boxcox

    Returns
    -------
    alternative distribution dataframes with either line or specimen as index
    """

    # Group by line and sequntaily run
    line_groupby = input_data.groupby('line')

    label_names = list(input_data.drop(['staging', 'line'], axis='columns').columns)

    baseline = input_data[input_data['line'] == 'baseline']
    baseline['genotype'] = 'wt'

    alt_line_pvalues = []
    alt_spec_pvalues = []

    # Get line-level alternative distributions
    for line_id, line_df in line_groupby:

        if line_id == 'baseline':
            continue

        line_df['genotype'] = 'hom'
        line_df.drop(['line'], axis=1)  # ?

        df_wt_mut = pd.concat([baseline, line_df])

        # Lm code needs daatpoints in numpy array and genotype+staging in dataframe
        data = df_wt_mut.drop(columns=['staging', 'line', 'genotype']).values
        info = df_wt_mut[['staging', 'line', 'genotype']]

        p, t = lm_r(data, info)  # returns p_values for all organs, 1 iteration
        res = [line_id] + list(p)
        alt_line_pvalues.append(res)

    # Get specimen-level alternative distributions
    mutants = input_data[input_data['line'] != 'baseline']
    # baselines = input_data[input_data['line'] == 'baseline']
    for specimen_id, row in mutants.iterrows():
        row['genotype'] = 'hom'
        line_id = row['line']
        df_wt_mut = baseline.append(row)
        data = df_wt_mut.drop(columns=['line', 'genotype', 'staging']).values
        info = df_wt_mut[['genotype', 'staging']]

        p, t = lm_r(data, info)  # returns p_values for all organs, 1 iteration
        res = [line_id, specimen_id] + list(p)
        alt_spec_pvalues.append(res)

    # result dataframes have either line or specimen in index then labels
    alt_line_df = pd.DataFrame.from_records(alt_line_pvalues, columns=['line'] + label_names, index='line')
    alt_spec_df = pd.DataFrame.from_records(alt_spec_pvalues, columns=['line', 'specimen'] + label_names, index='specimen')

    strip_x([alt_line_df, alt_spec_df])

    return alt_line_df, alt_spec_df
