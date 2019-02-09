#!/usr/bin/env python

"""
Given a Null and Altertiove distribution of a parameter (organ)
return a p-value threshold so that the false discovery will be set to 5%.
"""

import pandas as pd
import numpy as np

TARGET = 0.05


def get_thresholds(null_dist: pd.DataFrame, alt_dist: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate the per-organ p-value thresholds
    Given a wildtype null distribution of p-values and a alterntive (mutant)  distribution
    find the larget pvalue threshold that would give a FDR < 0.05

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
    dataframe
        columns: organs
        rows: p-value threshold

    """
    results = []

    for label in null_dist:

        wt_pvals = np.sort(null_dist[label].values)
        mut_pvals = np.sort(alt_dist[label].values)

        # TODO: the combined p-values should be sorted

        # Merge the p-values together get a list of available thresholds to use
        all_p = list(wt_pvals) + list(mut_pvals)
        all_p.sort()

        pthresh_fdrs = []

        # For every avaiable p-value from the null + alternative distributions, get the associated FDR
        for p_to_test in all_p:

            fdr_at_thresh = fdr_calc(wt_pvals, mut_pvals, p_to_test)

            if fdr_at_thresh:
                pthresh_fdrs.append((p_to_test, fdr_at_thresh))

        if pthresh_fdrs:

            # If we do not have a p < 0.05. So go for min FDR value
            if min(pthresh_fdrs, key=lambda x: x[0])[1] > 0.05:
                best_p, best_fdr = min(pthresh_fdrs, key=lambda x: x[1])  # x[1] is fdr

            # We have p values < 0.05 so choose the largest
            else:
                pthresh_fdrs = [x for x in pthresh_fdrs if x[1] <= 0.05]
                best_p, best_fdr = max(pthresh_fdrs, key=lambda x: x[0])  # x[0] is p

            # print(best_p, best_fdr)
            num_hits = len(mut_pvals[mut_pvals <= best_p])

        else:
            best_fdr = 1
            best_p = 'NA'
            num_hits = 0

        # TODO: what about if the labels are not numbers
        results.append([int(label), best_p, best_fdr, num_hits])

    header = ['label', 'p_thresh', 'fdr', 'num_hits_across_all_lines/specimens']
    result_df = pd.DataFrame.from_records(results, columns=header, index='label')
    result_df.sort_values(by='label', inplace=True)

    return result_df


def fdr_calc(wt_pvals, mut_pvals, thresh):
    ratio_wt_under_thresh = len(wt_pvals[wt_pvals < thresh]) / len(wt_pvals)

    ratio_mut_under_threshold = len(mut_pvals[mut_pvals < thresh]) / len(mut_pvals)

    try:
        fdr = ratio_wt_under_thresh / ratio_mut_under_threshold
    except ZeroDivisionError:
        # No mutants at this threshold. return a negative FDR
        return None
    return fdr






