#! /usr/bin/env python
import numpy as np
import scipy.stats.stats as stats
import sys
from os.path import join
import os
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common
import SimpleITK as sitk


def _flatten(arrays):
    one_d = []
    for arr in arrays:
        f = arr.flatten()  # try ravel to get a view rather than copy
        one_d.append(f)
    stacked = np.vstack(one_d)
    return stacked

def _get_data(paths):
    result = []
    for path in paths:
        img = sitk.ReadImage(path)
        array = sitk.GetArrayFromImage(img)
        array = array.astype(np.uint8)
        print path
        print array.shape
        result.append(array)
    return result

wt_dir = sys.argv[1]
mut_dir = sys.argv[2]

wt_paths = common.GetFilePaths(wt_dir)
mut_paths = common.GetFilePaths(mut_dir)

wt_data = _get_data(wt_paths)
mut_data = _get_data(mut_paths)

flat_mut = _flatten(mut_data)
flat_wt = _flatten(wt_data)

tstats, pvalues = stats.ttest_ind(flat_wt, flat_mut)

