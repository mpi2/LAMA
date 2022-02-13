#!/usr/bin/env python

"""
Given a Null and Altertiove distribution of a parameter (organ)
return a p-value threshold so that the false discovery will be set to 5%.
"""

import pandas as pd
import numpy as np

TESTING = False  # If set to true the p-threshold will be set high an dthe fdr < 0.05


# to get some positive hits for testing


def get_thresholds(null_dist: pd.DataFrame, alt_dist: pd.DataFrame, target_threshold: float = 0.05,
                   two_way: bool = False) -> pd.DataFrame:
    """
    Calculate the per-organ p-value thresholds
    Given a wild type null distribution of p-values and a alternative (mutant)  distribution
    find the largest p-value threshold that would give a FDR < 0.05

    Parameters
    ----------
    null_dist
        columns: organs
        rows: an iteration result

    alt_dist
        columns: organs
        rows: p-values from a mutant line or specimen

    target_threshold
        The target FDR threshold
    Returns
    -------
    pd.DataFrame
    Contains the pvalue threshold and associated fdr value for  each label
    Also has columns with value sthat were used to generate these values (added for debugging)
        columns:
            ['label', 'p_thresh', 'fdr', 'num_hits_across_all_lines/specimens',
            num_null, num_null_<=_thresh, num_alt, num_alt_<=_thresh]

    """
    results = []

    # Prevent overwriting originals
    null_dist = null_dist.copy()
    alt_dist = alt_dist.copy()

    for label in null_dist:

        if two_way:
            # TODO seem if you can improve performance
            # convert this back to an array

            wt_pvals = np.vstack(null_dist[label].values).transpose()
            mut_pvals = np.vstack(alt_dist[label].values).transpose()

            # Join array and vstack them - transpose so the effects are rows
            all_p = np.concatenate((wt_pvals, mut_pvals), axis=1)

            # sort for each effect

            for row in all_p:
                row.sort()

            # crete empty lists

            p_fdr_df = []

            # iterate for each effect
            for i, row in enumerate(all_p):
                pthresh_fdrs = []
                row = [x for x in row if x <= 0.05]
                # what to do if the rows are empty
                for p_to_test in row:
                    # basically index only compares the correct effect

                    fdr_at_thresh = fdr_calc(wt_pvals[i], mut_pvals[i], p_to_test)

                    if fdr_at_thresh is not None:
                        pthresh_fdrs.append((p_to_test, fdr_at_thresh))

                p_fdr = pd.DataFrame.from_records(pthresh_fdrs, columns=['p', 'fdr'])

                p_fdr_df.append(p_fdr)



            # enumerate for performance
            for i, p_fdr in enumerate(p_fdr_df):

                if len(p_fdr) > 0:
                    p_under_target_fdr = p_fdr[p_fdr.fdr <= target_threshold]

                    if len(p_under_target_fdr) < 1:
                        lowest_fdr_row = p_fdr.loc[p_fdr['fdr'].idxmin()]
                        p_thresh = lowest_fdr_row['p']
                        best_fdr = lowest_fdr_row['fdr']
                    else:
                        row = p_fdr.loc[p_under_target_fdr.p.idxmax()]
                        p_thresh = row['p']
                        best_fdr = row['fdr']

                    num_hits = len(mut_pvals[i][mut_pvals[i] <= p_thresh])
                    num_null = len(wt_pvals[i])
                    num_alt = len(mut_pvals[i])

                    num_null_lt_thresh = len(wt_pvals[i][wt_pvals[i] <= p_thresh])

                else:
                    best_fdr = 1
                    p_thresh = np.NAN
                    num_hits = 0
                    num_null, num_null_lt_thresh, num_alt = ['NA'] * 3

                effect_list = ['genotype', 'treatment', 'interaction']
                results.append([int(label), effect_list[i], p_thresh, best_fdr,
                                num_null, num_null_lt_thresh, num_alt, num_hits])


        else:
            wt_pvals = null_dist[label].values
            mut_pvals = alt_dist[label].values

            wt_pvals.sort()
            mut_pvals.sort()
            all_p = list(wt_pvals) + list(mut_pvals)
            all_p.sort()
            pthresh_fdrs = []
            # For every available p-value from the null + alternative distributions, That is lower than 0.05
            # get the associated FDR for that threshold
            all_p = [x for x in all_p if x <= 0.05]
            for p_to_test in all_p:

                fdr_at_thresh = fdr_calc(wt_pvals, mut_pvals, p_to_test)

                if fdr_at_thresh is not None:
                    pthresh_fdrs.append((p_to_test, fdr_at_thresh))

            # create a dataframe of p-value thresholds and associated fdrs
            p_fdr_df = pd.dataframe.from_records(pthresh_fdrs, columns=['p', 'fdr'])

            if len(p_fdr_df) > 0:
                p_under_target_fdr = p_fdr_df[p_fdr_df.fdr <= target_threshold]

                if len(p_under_target_fdr) < 1:
                    lowest_fdr_row = p_fdr_df.loc[p_fdr_df['fdr'].idxmin()]
                    p_thresh = lowest_fdr_row['p']
                    best_fdr = lowest_fdr_row['fdr']
                else:
                    row = p_fdr_df.loc[p_under_target_fdr.p.idxmax()]
                    p_thresh = row['p']
                    best_fdr = row['fdr']

                # Total number of parameters across all lines that are below our p-value threshold
                num_hits = len(mut_pvals[mut_pvals <= p_thresh])

                num_null = len(wt_pvals)
                num_alt = len(mut_pvals)

                num_null_lt_thresh = len(wt_pvals[wt_pvals <= p_thresh])
                # num_alt_lt_thresh = len(mut_pvals[mut_pvals <= p_thresh])

            else:
                best_fdr = 1
                p_thresh = np.NAN
                num_hits = 0
                num_null, num_null_lt_thresh, num_alt = ['NA'] * 3

            results.append([int(label), p_thresh, best_fdr,
                            num_null, num_null_lt_thresh, num_alt, num_hits])

            # iteration over each group needs to be done separately

            #     if n_accept_pvals < 1:
            #         # No acceptable p-value threshold for this label. Choose minimum fdr.
            #         lowest_fdr_row = [p_fdr.loc[p_fdr['fdr'].idxmin()] for p_fdr in p_fdr_df]
            #         p_thresh = [row['p'] for row in lowest_fdr_row]
            #         best_fdr = [row['fdr'] for row in lowest_fdr_row]
            #
            # elif two_way:
            #     # this is amazing if it works
            #     rows = []
            #     for i, p_fdr in enumerate(p_fdr_df):
            #         rows.append(p_fdr.loc[p_under_target_fdr[i].p.idxmax()])
            #     p_thresh = [row['p'] for row in rows]
            #     best_fdr = [row['fdr'] for row in rows]

        # TODO: what about if the labels are not numbers

    if two_way:
        header = ['label', 'effect', 'p_thresh', 'fdr',
                  'num_null', 'num_null_lt_thresh', 'num_alt', 'num_alt_lt_thresh']

        result_df = pd.DataFrame.from_records(results, columns=header, index='label')
        result_df.sort_values(by=['label','effect'], inplace=True)

    else:
        header = ['label', 'p_thresh', 'fdr',
                  'num_null', 'num_null_lt_thresh', 'num_alt', 'num_alt_lt_thresh']

        result_df = pd.DataFrame.from_records(results, columns=header, index='label')
        result_df.sort_values(by='label', inplace=True)



    return result_df


def fdr_calc(null_pvals, alt_pvals, thresh, two_way=False) -> float:
    """
    Calculate the False Discovery Rate for a given p-value threshold and a null and alternative distribution
    Parameters

    ----------
    wt_pvals
    mut_pvals
    thresh

    Returns
    -------
    fdr [0.0,1.0]
    or None if both ratio_wt_under_thresh / ratio_mut_under_threshold are zero (Currently looking into ths)

    """
    null_pvals = np.sort(null_pvals)

    alt_pvals = np.sort(alt_pvals)

    ratio_wt_under_thresh = len(null_pvals[null_pvals < thresh]) / len(null_pvals)

    ratio_mut_under_threshold = len(alt_pvals[alt_pvals < thresh]) / len(alt_pvals)



    try:
        # if there are no wild-types, the FDR is 0?
        fdr = ratio_wt_under_thresh / ratio_mut_under_threshold
    except ZeroDivisionError:
        # No mutants at this threshold.
        return None
    # If the null is skewed to the right, we might get FDR values greater than 1, which does not make sense

    fdr = np.clip(fdr, 0, 1)
    return fdr
