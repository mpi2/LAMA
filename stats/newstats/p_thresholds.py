#!/usr/bin/env python

"""
We have pvalue threshold per organ for the organ volume analysis. See lab book 060818
Now I would like to optimize the threshold rather than choosing from a small list of potential thresholds.

see lab book 250818
"""

import pandas as pd
import numpy as np

import scipy.optimize as optimize
from scipy.optimize.nonlin import NoConvergence
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

    # wt_43 = wt_df['43'].sort_values()
    # print(wt_43)  OK!

    for label in null_dist:

        wt_pvals = np.sort(null_dist[label].values)
        mut_pvals = np.sort(alt_dist[label].values)

        # Merge the p-values together get a list of available thresholds to use
        all_p = list(wt_pvals) + list(mut_pvals)

        label_results = []

        for p_to_test in all_p:
            # if label == '43':
            #
            #     print(wt_pvals)  # OK! up to here
            #     print(mut_pvals)

            fdr = fdr_calc(wt_pvals, mut_pvals, p_to_test)

            if fdr:
                label_results.append((p_to_test, fdr))

        # If we do not have a p < 0.05. So go for min FDR value
        if min(label_results, key=lambda x: x[0])[1] > 0.05:
            best_p, best_fdr = min(label_results, key=lambda x: x[1])  # x[1] is fdr

        # We have p values < 0.05 so choose the largest
        else:
            label_results = [x for x in label_results if x[1] <= 0.05]
            best_p, best_fdr = max(label_results, key=lambda x: x[0])  # x[0] is p

        # print(best_p, best_fdr)
        num_hits = len(mut_pvals[mut_pvals <= best_p])
        results.append([int(label), best_p, best_fdr, num_hits])

    header = ['label', 'p_thresh', 'fdr', 'num_hits_across_all_lines/specimens']
    result_df = pd.DataFrame.from_records(results, columns=header)
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






