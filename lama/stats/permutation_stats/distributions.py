#!/usr/bin/env python3

"""
06/08/18

* Determine a null distributiom in order to set an appropriate p-value threshold


Steps
-----

For each mutant line
    get n number of homs
    relabel n wildtypes as mutants
    run organ volume LM   organ_volume ~ genotype + 'staging metric'



Notes
-----
The p and t-values that are returned by lm() which calls R lm() usually has specimen-level values appended to the
matrix.

This does not happen here as mutants must be labelled 'mutant' for this to happen (see Rscrips/lmFast.R)
TDOD: If the relabelled baselines were relabelled as 'mutant' instead og 'hom' or 'synthetic_hom' the speciemn-level
p-values, tvalues could be obtained from the same call the lm() as in the standard stats module.

"""

from os.path import expanduser
from typing import Union, Tuple, List
import random
from pathlib import Path
import math
from collections import Counter

import pandas as pd
import numpy as np
from scipy.special import comb
import statsmodels.formula.api as smf
from joblib import Parallel, delayed
import datetime
from logzero import logger
from tqdm import tqdm
import random
import itertools

from lama.stats.linear_model import lm_r, lm_sm

home = expanduser('~')


def generate_random_combinations(data: pd.DataFrame, num_perms):
    logger.info('generating permutations')
    data = data.drop(columns='staging', errors='ignore')
    line_specimen_counts = get_line_specimen_counts(data)

    result = {}
    # now for each label calcualte number of combinations we need for each
    for label in line_specimen_counts:
        label_indices_result = []
        counts = line_specimen_counts[label].value_counts()

        number_of_lines = counts[counts.index != 0].sum() # Drop the lines with zero labels (have been qc'd out)
        ratios = counts[counts.index != 0] / number_of_lines
        num_combs = num_perms * ratios

        # get wt data for label
        label_data = data[[label, 'line']]
        label_data = label_data[label_data.line == 'baseline']
        label_data = label_data[~label_data[label].isna()]

        # Soret out the numnbers

        # Sort out the numbers
        records = []

        for n, n_combs_to_try in num_combs.items():
            n_combs_to_try = math.ceil(n_combs_to_try)
            max_combs = int(comb(len(label_data), n))
            # logger.info(f'max_combinations for n={n} and wt_n={len(label_data)} = {max_combs}')
            records.append([n, n_combs_to_try, max_combs])
        df = pd.DataFrame.from_records(records, columns=['n', 'num_combs', 'max_combs'], index='n').sort_index(ascending=True)

        # test whether it's possible to have this number of permutations with data structure
        print(f'Max combinations for label {label} is {df.max_combs.sum()}')
        if num_perms > df.max_combs.sum():
            raise ValueError(f'Max number of combinations is {df.max_combs.sum()}, you requested {num_perms}')

        # Now spread the overflow from any ns to other groups
        while True:
            df['overflow'] = df.num_combs - df.max_combs  # Get the 'overflow' num permutations over maximumum unique
            groups_full = df[df.overflow >= 0].index

            df['overflow'][df['overflow'] < 0] = 0
            extra = df[df.overflow > 0].overflow.sum()

            df.num_combs -= df.overflow

            if extra < 1: # All combimation amounts have been distributed
                break

            num_non_full_groups = len(df[df.overflow >= 0])

            top_up_per_group = math.ceil(extra / num_non_full_groups)
            for n, row in df.iterrows():
                if n in groups_full:
                    continue
                # Add the topup amount
                row.num_combs += top_up_per_group

        # now generate the indices
        indx = label_data.index
        for n, row in df.iterrows():

            combs_gen = itertools.combinations(indx, n)
            for i, combresult in enumerate(combs_gen):
                if i == row.num_combs:
                    break
                label_indices_result.append(combresult)

        result[label] = label_indices_result
    return result


def max_combinations(num_wts: int, line_specimen_counts: dict) -> int:
    """

    num_wts
        Total number of wild types
    lines_n
        {line_n: number of lines}
    Returns
    -------
    Maximum number or permutations
    """
    # calcualte the maximum number of permutations allowed given the WT n and the mutant n line structure.
    results = {}

    counts = line_specimen_counts.iloc[:,0].value_counts()
    for n, num_lines in counts.items():
        # Total number of combinations given WT n and this mut n
        total_combs_for_n = int(comb(num_wts, n))
        # Now weight based on how many lines have this n
        total_combs_for_n /= num_lines
        results[n] = total_combs_for_n

    return int(min(results.values()))


def get_line_specimen_counts(input_data: pd.DataFrame) -> pd.DataFrame:
    """
    For each mutant line get the number of specimens per label. Does not inlude specimen labels that are NAN thereby
    accounting for QC'd out labels

    Parameters
    ----------
    input_data
            index: specimen_id
            cols: label numbers (sppended with x for use with statsmodels) eg. x1, x2
                  line id e.g baseline

    Returns
    -------
    index: line_id
    cols: label numbers (sppended with x for use with statsmodels) eg. x1, x2


    """
    if 'line' in input_data:
        col = 'line'
    elif 'genotype' in input_data:
        col = 'genotype'
    line_specimen_counts = input_data[input_data[col] != 'baseline'].groupby(col).count()
    return line_specimen_counts


def null(input_data: pd.DataFrame,
         num_perm: int,) -> Tuple[pd.DataFrame, pd.DataFrame, List]:
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

    # Store p-value and t-value results. One tuple (len==num labels) per iteration
    spec_p = []

    # Create synthetic specimens by iteratively relabelling each baseline as synthetic mutant
    baselines = input_data[input_data['line'] == 'baseline']

    # Get the line specimen n numbers. Keep the first column
    # line_specimen_counts = get_line_specimen_counts(input_data)
    # Pregenerate all the combinations
    wt_indx_combinations = generate_random_combinations(input_data, num_perm)

    # Split data into a numpy array of raw data and dataframe for staging and genotype fpr the LM code
    data = baselines.drop(columns=['staging', 'line']).values
    info = baselines[['staging', 'line']]

    # Get the specimen-level null distribution. i.e. the distributuion of p-values obtained from relabelling each
    # baseline once. Loop over each specimen and set to 'synth_hom'
    for index, _ in info.iterrows():
        info.loc[:, 'genotype'] = 'wt'               # Set all genotypes to WT
        info.loc[[index], 'genotype'] = 'synth_hom'  # Set the ith baseline to synth hom
        row = data[info.index.get_loc(index), :]

        # Get columns (labels) where the mutant specimen (as it's line level)
        # has a null value (i.e. This specimen is QC-flagged at these labels)
        # Set all values to zero
        labels_to_skip = np.isnan(row)
        if any(labels_to_skip):
            # Set the whole label column to zero. R:lm() will return NaN for this column
            d = np.copy(data)
            d[:, labels_to_skip] = 0.0
        else:
            d = data

        # Get a p-value for each organ
        p, t = lm_sm(d, info)  # TODO: move this to statsmodels

        # Check that there are equal amounts of p-values than there are data points
        if len(p) != data.shape[1]:
            raise ValueError(f'The length of p-values results: {data.shape[1]} does not match the length of the input data: {len(p)}')

        spec_p.append(p)

    spec_df = pd.DataFrame.from_records(spec_p, columns=label_names)

    line_df = null_line(wt_indx_combinations, baselines, num_perm)

    return strip_x([line_df, spec_df])


def null_line(wt_indx_combinations: dict,
              data: pd.DataFrame,
              num_perms=1000) -> pd.DataFrame:
    """
    Generate pvalue null distributions for all labels in 'data'
    NaN values are excluded potentailly resultnig in different sets of specimens for each label. This makes it tricky to
    use Rs vectorised lm() function, so we use statsmodels here and just loop over the data, with each label in its own
    process using joblib.

    Parameters
    ----------
    line_specimen_counts
        The number of specimens in each mutant line. Used to specify the number of synthetic mutants so the null
        matches our data in terms of sample number
    data
        Label data in each column except last 2 which are 'staging' and 'genotype'
    num_perms
        Usually about 10000

    Returns
    -------
    DataFrame of null distributions. Each label in a column

    Notes
    -----
    If QC has been applied to the data, we may have some NANs
    """

    def prepare(label):
        return data[[label, 'staging', 'genotype']]

    data = data.rename(columns={'line': 'genotype'})

    starttime = datetime.datetime.now()

    cols = list(data.drop(['staging', 'genotype'], axis='columns').columns)

    # Run each label on a thread
    pdists = Parallel(n_jobs=-1)(delayed(_null_line_thread)
                                 (prepare(i), num_perms, wt_indx_combinations, i) for i in tqdm(cols))

    line_pdsist_df = pd.DataFrame(pdists).T
    line_pdsist_df.columns = cols

    endtime = datetime.datetime.now()
    elapsed = endtime - starttime
    print(f'Time taken for null distribution calculation: {elapsed}')
    return line_pdsist_df


def _null_line_thread(*args) -> List[float]:
    """
    Create a null distribution for a single label.  This can put put onto a thread or process

    Returns
    -------
    pvalue distribution
    """
    data, num_perms, wt_indx_combinations, label = args
    print('Generating null for', label)

    label = data.columns[0]

    data = data.astype({label: np.float,
                        'staging': np.float})

    synthetics_sets_done = []

    line_p = []

    perms_done = 0

    # Get combinations of WT indices for current label
    indxs = wt_indx_combinations[label]
    for comb in indxs:
        data.loc[:, 'genotype'] = 'wt'
        data.loc[data.index.isin(comb), 'genotype'] = 'synth_hom'
        # _label_synthetic_mutants(data, n, synthetics_sets_done)

        perms_done += 1

        model = smf.ols(formula=f'{label} ~ C(genotype) + staging', data=data, missing='drop')
        fit = model.fit()
        p = fit.pvalues['C(genotype)[T.wt]']

        line_p.append(p)
    return line_p


def _label_synthetic_mutants(info: pd.DataFrame, n: int, sets_done: List) -> bool:
    """
    Given a dataframe of wild type data, relabel n baselines as synthetic mutant in place.
    Keep track of combinations done in sets_done and do not duplicate

    Parameters
    ----------
    info
        columns
            label_num with 'x' prefix e.g. 'x1'
            staging
            genotype
    n
        how many specimens to relabel
    sets_done
        Contains Sets of previously selected specimen IDs

    Returns
    -------
    True if able to find a new combination of mutant and wildtype
    False if no more combinations are available for that specific n
    """

    # Set all to wt genotype
    info.loc[:, 'genotype'] = 'wt'

    # label n number of baselines as mutants from the maximum number of combination

    # max_comb = int(comb(len(info), n))
    #
    # for i in range(max_comb):
    #     synthetics_mut_indices = random.sample(range(0, len(info)), n)
    #     i += 1
    #     if not set(synthetics_mut_indices) in sets_done:
    #         break
    #
    # if i > max_comb - 1:
    #     msg = f"""Cannot find unique combinations of wild type baselines to relabel as synthetic mutants
    #     With a baseline n  of {len(info)}\n. Choosing {n} synthetics.
    #     Try increasing the number of baselines or reducing the number of permutations"""
    #     logger.warn(msg)
    #     raise ValueError(msg)

    sets_done.append(set(synthetics_mut_indices))

    # info.ix[synthetics_mut_indices, 'genotype'] = 'synth_hom'  # Why does iloc not work here?
    info.loc[info.index[synthetics_mut_indices], 'genotype'] = 'synth_hom'
    return True


def strip_x(dfs):
    for df in dfs:
        df.columns = [x.strip('x') for x in df.columns]
        yield df


def alternative(input_data: pd.DataFrame,
                plot_dir: Union[None, Path] = None,
                boxcox: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Generate alterntive (mutant) distributions for line and pecimen-level data

    Parameters
    ----------
    input_data
    plot_dir
    boxcox

    Returns
    -------
    alternative distribution dataframes with either line or specimen as index.
        0: line-level p values
        1: specimen-level p values
        2: line-level t-values
        3: specimen-level t-values
    """

    # Group by line and sequntaily run
    info_columns = ['staging', 'line', 'genotype']  # Columns of non-organ volumes in input_data

    line_groupby = input_data.groupby('line')

    label_names = list(input_data.drop(['staging', 'line'], axis='columns').columns)

    baseline = input_data[input_data['line'] == 'baseline']
    baseline.loc[:, 'genotype'] = 'wt'

    alt_line_pvalues = []
    alt_spec_pvalues = []
    alt_line_t = []
    alt_spec_t = []

    # Get line-level alternative distributions
    for line_id, line_df in line_groupby:

        if line_id == 'baseline':
            continue

        line_df.loc[:, 'genotype'] = 'hom'
        line_df.drop(['line'], axis=1)

        df_wt_mut = pd.concat([baseline, line_df])

        # Lm code needs daatpoints in numpy array and genotype+staging in dataframe
        data_df = df_wt_mut.drop(columns=info_columns)

        # Get columns (labels) where all speciemns have null values (i.e. all the line is QC-flagged at these labels)
        labels_to_skip = [col for col, isany in line_df.any().iteritems() if not isany]
        if labels_to_skip:
            # Set the whole label column to zero. R:lm() will return NaN for this column
            data_df[labels_to_skip] = 0.0

        # Get a numpy array of the organ volumes
        data = data_df.values

        info = df_wt_mut[['staging', 'line', 'genotype']]

        p: np.array
        t: np.array

        p, t = lm_sm(data, info)  # returns p_values for all organs, 1 iteration

        res_p = [line_id] + list(p)  # line_name, label_1, label_2 ......
        alt_line_pvalues.append(res_p)

        res_t = [line_id] + list(t)
        alt_line_t.append(res_t)

    ### Get specimen-level alternative distributions ###
    mutants = input_data[input_data['line'] != 'baseline']
    # baselines = input_data[input_data['line'] == 'baseline']
    for specimen_id, row in mutants.iterrows():
        row['genotype'] = 'hom'
        line_id = row['line']
        df_wt_mut = baseline.append(row)
        data_df = df_wt_mut.drop(columns=['line', 'genotype', 'staging'])

        # Get columns (labels) where the mutant specimen (as it's line level)
        # has a null value (i.e. This specimen is QC-flagged at these labels)
        labels_to_skip = [col for col, isany in pd.DataFrame(row).T.any().iteritems() if not isany]
        if labels_to_skip:
            # Set the whole label column to zero. R:lm() will return NaN for this column
            data_df[labels_to_skip] = 0.0

        data = data_df.values

        info = df_wt_mut[['genotype', 'staging']]

        p, t = lm_sm(data, info)  # returns p_values for all organs, 1 iteration
        res_p = [line_id, specimen_id] + list(p)
        alt_spec_pvalues.append(res_p)

        res_t = [specimen_id] + list(t)
        alt_spec_t.append(res_t)

    # result dataframes have either line or specimen in index then labels
    alt_line_df = pd.DataFrame.from_records(alt_line_pvalues, columns=['line'] + label_names, index='line')
    alt_spec_df = pd.DataFrame.from_records(alt_spec_pvalues, columns=['line', 'specimen'] + label_names,
                                            index='specimen')

    alt_line_t_df = pd.DataFrame.from_records(alt_line_t, columns=['line'] + label_names, index='line')
    alt_spec_t_df = pd.DataFrame.from_records(alt_spec_t, columns=['specimen'] + label_names,
                                            index='specimen')

    return strip_x([alt_line_df, alt_spec_df, alt_line_t_df, alt_spec_t_df])
