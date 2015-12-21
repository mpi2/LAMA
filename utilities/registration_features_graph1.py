#! /usr/bin/env python

"""

"""

import SimpleITK as sitk
from SimpleITK import ReadImage, WriteImage, GetArrayFromImage, GetImageFromArray
import numpy as np
from os.path import join
import sys
import os
import matplotlib.pyplot as plt
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common


def run(wt_dir, mut_dir, title):


    for i, wt_path in enumerate(common.GetFilePaths(wt_dir)):
        # Try out gettiong features along z axis

        img = ReadImage(wt_path)
        arr = GetArrayFromImage(img)
        if i == 0:
            wt_mean = arr
        else:
            wt_mean += arr

    for i, mut_path in enumerate(common.GetFilePaths(mut_dir)):
        # Try out gettiong features along z axis

        img = ReadImage(mut_path)
        arr = GetArrayFromImage(img)
        if i == 0:
            mut_mean = arr
        else:
            mut_mean += arr



    plt.plot(list(np.mean(wt_mean, axis=(1,2))), label="Control")
    plt.plot(list(np.mean(mut_mean, axis=(1,2))), label="Treated")
    plt.title(title)
    plt.legend(loc='upper left')
    plt.show()

if __name__ == '__main__':
    wt_dir = sys.argv[1]
    mut_dir = sys.argv[2]
    title = sys.argv[3]
    run(wt_dir, mut_dir, title)
