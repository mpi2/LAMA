#! /usr/bin/env python
import numpy as np
import scipy.stats.stats as stats
import sys
from os.path import join
import os
import SimpleITK as sitk


def GetFilePaths(folder, extension_tuple=('.nrrd', '.tiff', '.tif', '.nii', '.bmp', 'jpg', 'mnc', 'vtk'), pattern=None):
    """
    Test whether input is a folder or a file. If a file or list, return it.
    If a dir, return all images within that directory.
    Optionally test for a pattern to sarch for in the filenames
    """
    if not os.path.isdir(folder):
        if isinstance(folder, basestring):
            return [folder]
        else:
            return folder
    else:
        paths = []
        for root, _, files in os.walk(folder):
            for filename in files:
                if filename.lower().endswith(extension_tuple):
                    if pattern:
                        if pattern and pattern not in filename:
                            continue
                    #paths.append(os.path.abspath(os.path.join(root, filename))) #Broken on shared drive
                    paths.append(os.path.join(root, filename))
        return paths

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

wt_paths = GetFilePaths(wt_dir)
mut_paths = GetFilePaths(mut_dir)

wt_data = _get_data(wt_paths)
mut_data = _get_data(mut_paths)

flat_mut = _flatten(mut_data)
flat_wt = _flatten(wt_data)

tstats, pvalues = stats.ttest_ind(flat_wt, flat_mut)

