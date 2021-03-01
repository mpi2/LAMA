import random

import pandas as pd
import numpy as np

from lama.stats.permutation_stats import distributions
random.seed(1)

def test_line_null():
    """
    Test the code to generate null distributions at the line level
    """
    # data = pd.DataFrame.from_dict({
    #     'x1_label': np.random.randn(100),
    #     'x2_label': np.random.randn(100),
    #     'staging': np.random.randn(100),
    #     'genotype': ['wt'] * 100
    # }, orient='index').T
    #
    num_labels = 2
    num_wt_specimens = 100
    num_perms = 100

    d = {}
    for i in range(num_labels):
        l = f'x{i}'
        d[l] = np.random.randn(num_wt_specimens)
    d['staging'] = np.random.randn(num_wt_specimens)
    d['line'] = ['baseline'] * num_wt_specimens
    data = pd.DataFrame.from_dict(d, orient='index').T



    # Set one of the speciemns' first label to NaN
    data.iloc[4, 0] = pd.NA

    specime_line_nums = np.random.randint(1, 8, 20)  # Twenty lines with a mutant sample size between 1 and 8
    null_dist = distributions.null_line(specime_line_nums,
                                        data,
                                        num_perms)