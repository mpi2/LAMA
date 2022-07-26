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


def random_combination(iterable, r):
    "Random selection from itertools.combinations(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(range(n), r))
    return tuple(pool[i] for i in indices)

def recursive_comb_maker(lst, n, steps, i, recurs_results):
    """makes a combination from a list given the n_groups"""
    # generate first set of combinations
    #combs_gen = list(itertools.combinations(lst, n // (steps + 1)))
    # Randomly choose a combination
    #comb_result = random.choices(combs_gen, k=1)
    comb_result = random_combination(lst, n // (steps + 1))
    #comb_result = random.choices(combs_gen, k=1)
    # append result into a list
    recurs_results.append(comb_result)
    #print(list(lst))
    # break if you've done the right amount of steps or you have a lot of combinations, else call it recurvisely
    if i == steps:
        return recurs_results
    else:
        return recursive_comb_maker(list(filter(lambda val: val not in comb_result[0], lst)), n, steps, i + 1, recurs_results)


def generate_random_combinations(data: pd.DataFrame, num_perms):
    logger.info('generating permutations')
    data = data.drop(columns='staging', errors='ignore')
    line_specimen_counts = get_line_specimen_counts(data)

    result = {}
    # now for each label calcualte number of combinations we need for each
    for label in line_specimen_counts:
        label_indices_result = []
        counts = line_specimen_counts[label].value_counts()

        number_of_lines = counts[counts.index != 0].sum()  # Drop the lines with zero labels (have been qc'd out)
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
        df = pd.DataFrame.from_records(records, columns=['n', 'num_combs', 'max_combs'], index='n').sort_index(
            ascending=True)

        # test whether it's possible to have this number of permutations with data structure
        print(f'Max combinations for label {label} is {df.max_combs.sum()}')
        if num_perms > df.max_combs.sum():
            raise ValueError(f'Max number of combinations is {df.max_combs.sum()}, you requested {num_perms}')

        # Now spread the overflow from any ns to other groups

        # Kyle - What is this?
        while True:
            df['overflow'] = df.num_combs - df.max_combs  # Get the 'overflow' num permutations over maximumum unique
            groups_full = df[df.overflow >= 0].index

            df['overflow'][df['overflow'] < 0] = 0
            extra = df[df.overflow > 0].overflow.sum()

            df.num_combs -= df.overflow

            if extra < 1:  # All combimation amounts have been distributed
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


def generate_random_two_way_combinations(data: pd.DataFrame, num_perms):
    logger.info('generating permutations')
    data = data.drop(columns='staging', errors='ignore')
    line_specimen_counts = get_line_specimen_counts(data, two_way=True)
    n_groups = get_two_way_n_groups(data)

    result = {}
    # now for each label calculate number of combinations we need for each
    for label in line_specimen_counts:
        label_indices_result = []
        counts = line_specimen_counts[label].value_counts()

        number_of_lines = counts[counts.index != 0].sum()  # Drop the lines with zero labels (have been qc'd out)
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
            max_combs = two_way_max_combinations(n, n_groups)
            # logger.info(f'max_combinations for n={n} and wt_n={len(label_data)} = {max_combs}')
            records.append([n, n_combs_to_try, max_combs])
        df = pd.DataFrame.from_records(records, columns=['n', 'num_combs', 'max_combs'], index='n').sort_index(
            ascending=True)

        # test whether it's possible to have this number of permutations with data structure
        # so for the two-way this wont change
        #print(f'Max combinations for label {label} is {max_combs}')
        if num_perms > df.max_combs.sum():
            raise ValueError(f'Max number of combinations is {max_combs}, you requested {num_perms}')


        indx = label_data.index

        label_indices_result = []

        # generate combinations for no. of desired perms
        for perm in range(num_perms):
            full_combination = recursive_comb_maker(indx, len(indx), n_groups - 1, i=1, recurs_results=[])
            # print(full_combination)

            label_indices_result.append(full_combination)
            # reset result list
            recurs_results = []

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

    counts = line_specimen_counts.iloc[:, 0].value_counts()
    for n, num_lines in counts.items():
        # Total number of combinations given WT n and this mut n
        total_combs_for_n = int(comb(num_wts, n))
        # Now weight based on how many lines have this n
        total_combs_for_n /= num_lines
        results[n] = total_combs_for_n

    return int(min(results.values()))


def two_way_max_combinations(num_wts: int, n_groups: int) -> int:
    """So the number of max combinations two the two-way stuff is completely different
    basically, its PI (product notation) i = 0 -> 3 [n_wts-(n_wts/4)* i C n_wts/4]
    once again  - I don't think I care about lines right now"""
    # I'll leave the line architecture alive for now

    # Total number of combinations given WT n and this mut n
    # I dont need n - you can get the total per group from the wts
    n_per_group = num_wts // n_groups
    comb_per_group = [(comb(num_wts - (n_per_group * i), n_per_group)) for i in range(0, n_groups)]
    total_combs_for_n = np.product(comb_per_group)
    # Now weight based on how many lines have this n
    return int(total_combs_for_n)


def get_two_way_n_groups(input_data: pd.DataFrame) -> int:
    """just get the number of different groups for two_way"""
    return len(input_data.groupby('line').count())


def get_line_specimen_counts(input_data: pd.DataFrame, two_way: bool = False) -> pd.DataFrame:
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
    if two_way:
        # So for the two-way I need the wild-type counts for the null distribution
        line_specimen_counts = input_data[input_data[col] == 'baseline'].groupby(col).count()
    else:
        line_specimen_counts = input_data[input_data[col] != 'baseline'].groupby(col).count()
    return line_specimen_counts


def null(input_data: pd.DataFrame,
         num_perm: int, two_way: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame, List]:
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

    two_way
        makes it a two-way null

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
    wt_indx_combinations = generate_random_two_way_combinations(input_data, num_perm) if two_way \
        else generate_random_combinations(input_data, num_perm)

    # Split data into a numpy array of raw data and dataframe for staging and genotype fpr the LM code
    data = baselines.drop(columns=['staging', 'line']).values
    info = baselines[['staging', 'line']]

    # Get the specimen-level null distribution. i.e. the distributuion of p-values obtained from relabelling each
    # baseline once. Loop over each specimen and set to 'synth_hom'
    if two_way:  # two-way flag
        # so for two_way we have three spec tests - geno, treat and int
        geno_info = info.copy()
        treat_info = info.copy()
        inter_info = info.copy()
        for index, meta in info.iterrows():

            geno_info.loc[:, 'genotype'] = 'wt'  # Set all genotypes to WT
            geno_info.loc[:, 'treatment'] = 'veh'  # Set all treatments to vehicle

            treat_info.loc[:, 'genotype'] = 'wt'  # Set all genotypes to WT
            treat_info.loc[:, 'treatment'] = 'veh'  # Set all treatments to vehicle

            inter_info.loc[:, 'genotype'] = 'wt'  # Set all genotypes to WT
            inter_info.loc[:, 'treatment'] = 'veh'  # Set all treatments to vehicle
            

            # now relabel based on the combinations

            #for geno effect label index as just synth mut:
            geno_info.loc[[index], 'genotype'] = 'synth_mut'

            #similarly for treatment:
            treat_info.loc[[index], 'treatment'] = 'synth_treat'

            # So the first logical thing is synth index as the interaction
            inter_info.loc[[index], 'genotype'] = 'synth_mut'
            inter_info.loc[[index], 'treatment'] = 'synth_treat'
            # then randomly assign 1/3 as mut and 1/3 as treat
            # generate list but fuck the index off:
            all_rows = info.drop(index).index.values

            #should be a list of two
            combs = recursive_comb_maker(all_rows,n=len(all_rows),steps=2,i=1,recurs_results=[])

            #assign synth int mutants
            inter_info.loc[combs[0][0], 'genotype'] = 'synth_mut'

            #assign synth int treatments
            inter_info.loc[combs[1][0], 'treatment'] = 'synth_treat'

            for _info in [geno_info,treat_info,inter_info]:
                row = data[_info.index.get_loc(index), :]
                labels_to_skip = np.isnan(row)
                if any(labels_to_skip):
                    # Set the whole label column to zero. R:lm() will return NaN for this column
                    d = np.copy(data)
                    d[:, labels_to_skip] = 0.0
                else:
                    d = data
                p,t = lm_sm(d, _info, two_way=True)  # TODO: move this to statsmodels
                if len(p) != data.shape[1]:
                    raise ValueError(
                        f'The length of p-values results: {data.shape[1]} does not match the length of the input data: {len(p)}')

                spec_p.append(p)

    else:
        for index, _ in info.iterrows():
            info.loc[:, 'genotype'] = 'wt'  # Set all genotypes to WT
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
                raise ValueError(
                    f'The length of p-values results: {data.shape[1]} does not match the length of the input data: {len(p)}')

            spec_p.append(p)

    spec_df = pd.DataFrame.from_records(spec_p, columns=label_names)

    line_df = null_line(wt_indx_combinations, baselines, num_perm, two_way=two_way)
    #print(line_df['x3'])
    return strip_x([line_df, spec_df])


def null_line(wt_indx_combinations: dict,
              data: pd.DataFrame,
              num_perms=1000, two_way:bool=False) -> pd.DataFrame:
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
        if two_way:
            return data[[label, 'staging', 'genotype','treatment']]
        else:
            return data[[label, 'staging', 'genotype']]

    if two_way:
        #make genotype and treatment cols
        data['genotype'] = 'wt'
        data['treatment'] = 'veh'
        data = data.drop(['line'], axis='columns')
    else:
        data = data.rename(columns={'line': 'genotype'})

    starttime = datetime.datetime.now()

    cols = list(data.drop(['staging', 'genotype', 'treatment'], axis='columns').columns) if two_way else \
        list(data.drop(['staging', 'genotype'], axis='columns').columns)

    # Run each label on a thread
    if two_way:
        pdists = Parallel(n_jobs=-1)(delayed(_two_way_null_line_thread)
                                     (prepare(i), num_perms, wt_indx_combinations, i) for i in tqdm(cols))
    else:
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


def _two_way_null_line_thread(*args) -> List[float]:
    """
    same as _null_line_thread but for two way
    TODO: merge two_way and null_threads
    Parameters
    ----------
    args

    Returns
    -------

    """
    data, num_perms, wt_indx_combinations, label = args
    # print('Generating null for', label)

    label = data.columns[0]

    data = data.astype({label: np.float,
                        'staging': np.float})

    synthetics_sets_done = []

    line_p = []

    perms_done = 0

    # Get combinations of WT indices for current label
    indxs = wt_indx_combinations[label]
    for comb in indxs:
        # set up genotype and treatment
        data.loc[:, 'genotype'] = 'wt'
        data.loc[:,'treatment'] = 'veh'

        # mains
        data.loc[data.index.isin(comb[0]), 'genotype'] = 'synth_mut'
        data.loc[data.index.isin(comb[1]), 'treatment'] = 'synth_treat'

        # interactions
        data.loc[data.index.isin(comb[2]), 'genotype'] = 'synth_mut'
        data.loc[data.index.isin(comb[2]), 'treatment'] = 'synth_treat'

        # _label_synthetic_mutants(data, n, synthetics_sets_done)

        perms_done += 1

        fit = smf.ols(formula=f'{label} ~ genotype * treatment + staging', data=data, missing='drop').fit()
        # get all pvals except intercept and staging

        # fit.pvalues is a series - theefore you have to use .index
        p = fit.pvalues[~fit.pvalues.index.isin(['Intercept','staging'])]
        #pvalues go in the order of genotype, treatment, interaction.
        line_p.append(p.values)
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
                boxcox: bool = False, two_way: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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
    info_columns =['staging', 'line', 'genotype', 'treatment'] if two_way \
        else ['staging', 'line', 'genotype']  # Columns of non-organ volumes in input_data

    line_groupby = input_data.groupby('line')



    label_names = list(input_data.drop(['staging', 'line'], axis='columns').columns)

    baseline = input_data[input_data['line'] == 'baseline']
    baseline.loc[:, 'genotype'] = 'wt'

    alt_line_pvalues = []
    alt_spec_pvalues = []
    alt_line_t = []
    alt_spec_t = []

    # Get line-level alternative distributions
    two_way_df = pd.DataFrame()
    if two_way:
        # create the proper two-by-two groups
        for line_id, line_df in line_groupby:
            if line_id == 'baseline':
                line_df.loc[:, 'genotype'] = 'wt'
                line_df.loc[:, 'treatment'] = 'veh'
            elif line_id == 'mutants':
                line_df.loc[:, 'genotype'] = 'mut'
                line_df.loc[:, 'treatment'] = 'veh'
            if line_id == 'treatment':
                line_df.loc[:, 'genotype'] = 'wt'
                line_df.loc[:, 'treatment'] = 'treat'
            if line_id == 'mut_treat':
                line_df.loc[:, 'genotype'] = 'mut'
                line_df.loc[:, 'treatment'] = 'treat'
            # merge the labels
            two_way_df=pd.concat([two_way_df, line_df])

        two_way_df.drop(['line'], axis=1)
        data_df = two_way_df.drop(columns=info_columns)

        # Get columns (labels) where all speciemns have null values (i.e. all the line is QC-flagged at these labels)
        labels_to_skip = [col for col, isany in line_df.any().iteritems() if not isany]

        if labels_to_skip:
            # Set the whole label column to zero. R:lm() will return NaN for this column
            data_df[labels_to_skip] = 0.0

        data = data_df.values

        info = two_way_df[['staging', 'treatment', 'genotype']]

        p: np.array
        t: np.array

        p, t = lm_sm(data, info, two_way=True)

        res_p = ['two_way'] + list(p)  # line_name, label_1, label_2 ......
        alt_line_pvalues.append(res_p)

        res_t = ['two_way'] + list(t)
        alt_line_t.append(res_t)

        ### Great now the specimen values ###
        # split the df again? ##TODO make this less stupid?
        two_grped = two_way_df.groupby(['genotype','treatment'])

        baseline = two_grped.get_group(('wt','veh'))
        mutants = two_grped.get_group(('mut','veh'))
        treatment = two_grped.get_group(('wt','treat'))
        interaction = two_grped.get_group(('mut','treat'))

        def get_effects(group, inter=False):
            for specimen_id, row in group.iterrows():
                line_id = 'two_way'
                df_wt_mut = baseline.append(row)
                if inter:
                    df_wt_mut = df_wt_mut.append(mutants)
                    df_wt_mut = df_wt_mut.append(treatment)

                data_df = df_wt_mut.drop(columns=['line', 'genotype', 'treatment', 'staging'])

                # Get columns (labels) where the mutant specimen (as it's line level)
                # has a null value (i.e. This specimen is QC-flagged at these labels)
                labels_to_skip = [col for col, isany in pd.DataFrame(row).T.any().iteritems() if not isany]
                if labels_to_skip:
                    # Set the whole label column to zero. R:lm() will return NaN for this column
                    data_df[labels_to_skip] = 0.0

                data = data_df.values

                info = df_wt_mut[['genotype', 'treatment', 'staging']]
                #print("data", data)
                #print("info", info)
                #print("type of info", type(info))
                #print('Type of data', type(data[-1]), data[-1][0])
                p, t = lm_sm(data, info, two_way=True)  # returns p_values for all organs, 1 iteration

                res_p = [line_id, specimen_id] + list(p)
                alt_spec_pvalues.append(res_p)

                res_t = [specimen_id] + list(t)
                alt_spec_t.append(res_t)

        # genotype effect
        get_effects(mutants)

        # treatment effect
        get_effects(treatment)

        get_effects(interaction, inter=True)

    else:
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
