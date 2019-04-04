#!/usr/bin/env python

"""
Given a Null and Altertiove distribution of a parameter (organ)
return a p-value threshold so that the false discovery will be set to 5%.
"""

import pandas as pd
import numpy as np

TARGET = 0.05

TESTING = False  # If set to true the p-threshold will be set high an dthe fdr < 0.05
# to get some positive hits for testing


def get_thresholds(null_dist: pd.DataFrame, alt_dist: pd.DataFrame) -> pd.DataFrame:
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

    for label in null_dist:

        wt_pvals = null_dist[label].values
        mut_pvals = alt_dist[label].values

        # Merge the p-values together get a list of available thresholds to use
        all_p = list(wt_pvals) + list(mut_pvals)

        all_p.sort()

        pthresh_fdrs = []

        # For every available p-value from the null + alternative distributions,
        # get the associated FDR for that threshold
        for p_to_test in all_p:

            fdr_at_thresh = fdr_calc(wt_pvals, mut_pvals, p_to_test)

            if fdr_at_thresh:
                pthresh_fdrs.append((p_to_test, fdr_at_thresh))

        if pthresh_fdrs:

            # If we do not have a p < 0.05 go for min FDR value
            if min(pthresh_fdrs, key=lambda x: x[0])[1] > 0.05:
                best_p, best_fdr = min(pthresh_fdrs, key=lambda x: x[1])  # x[1] is fdr

            # We have p values < 0.05 so choose the largest p=value
            else:
                pthresh_fdrs = [x for x in pthresh_fdrs if x[1] <= 0.05]
                best_p, best_fdr = max(pthresh_fdrs, key=lambda x: x[0])  # x[0] is p

            num_hits = len(mut_pvals[mut_pvals <= best_p])

            num_null = len(wt_pvals)
            num_alt = len(mut_pvals)

            num_null_lt_thresh = len(wt_pvals[wt_pvals < best_p])
            num_alt_lt_thresh = len(mut_pvals[mut_pvals < best_p])

        else:
            best_fdr = 1
            best_p = np.NAN
            num_hits = 0
            num_null, num_null_lt_thresh, num_alt, num_alt_lt_thresh = ['NA'] * 4

        if TESTING:
            print('Warning testing pthresholds')
            best_p = 1.0
            best_fdr = 0.049

        # TODO: what about if the labels are not numbers
        results.append([int(label), best_p, best_fdr, num_hits,
                        num_null, num_null_lt_thresh, num_alt, num_alt_lt_thresh])

    header = ['label', 'p_thresh', 'fdr', 'num_hits_across_all_lines/specimens',
              'num_null', 'num_null_lt_thresh', 'num_alt', 'num_alt_lt_thresh']

    result_df = pd.DataFrame.from_records(results, columns=header, index='label')
    result_df.sort_values(by='label', inplace=True)

    return result_df


def fdr_calc(null_pvals, alt_pvals, thresh) -> float:
    """
    Calculate the False Discovery Rate for a given p-valuie threshold and a null and alternative distribution
    Parameters

    ----------
    wt_pvals
    mut_pvals
    thresh

    Returns
    -------
    fdr [0.0,1.0]

    """
    null_pvals = np.sort(null_pvals)
    alt_pvals = np.sort(alt_pvals)
    ratio_wt_under_thresh = len(null_pvals[null_pvals < thresh]) / len(null_pvals)

    ratio_mut_under_threshold = len(alt_pvals[alt_pvals < thresh]) / len(alt_pvals)

    try:
        fdr = ratio_wt_under_thresh / ratio_mut_under_threshold
    except ZeroDivisionError:
        # No mutants at this threshold.
        return None
    return fdr






