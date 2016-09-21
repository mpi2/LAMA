#!/usr/bin/env python

import os
import SimpleITK as sitk
import numpy as np
from matplotlib import pyplot as plt
from tempfile import TemporaryFile

EXTENSION = ('nrrd', 'tif')
XLIM = (-0.25, 0.15)
YLIM = (0, 700)


def run(wt_dir, mut_dir, out_dir):

    wt_paths = [os.path.join(wt_dir, x) for x in os.listdir(wt_dir) if x.endswith(EXTENSION)]
    mut_paths = [os.path.join(mut_dir, x) for x in os.listdir(mut_dir) if x.endswith(EXTENSION)]
    paths = wt_paths + mut_paths

    memmap_arrays = memmap_data(paths)

    mean = np.mean(memmap_arrays, axis=0)
    for i, array in enumerate(memmap_arrays):
        errors = mean - np.array(array)
        bins = np.linspace(XLIM[0], XLIM[1], num=2048)
        hist1, bins = np.histogram(errors, bins=bins)
        width = 1 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist1, width=width, color='blue', align='center', alpha=0.4, linewidth=0)
        title = os.path.splitext(os.path.basename(paths[i]))[0]
        plt.title(title)
        plt.xlim(XLIM)
        plt.ylim(YLIM)
        plt.savefig(os.path.join(out_dir, title + '.png'))
        plt.close()


def memmap_data(paths):
    mm_arrays = []
    for f in paths:
        img = sitk.ReadImage(f)
        arr = sitk.GetArrayFromImage(img)
        t = TemporaryFile()
        mm = np.memmap(t, dtype=arr.dtype, mode='w+', shape=arr.shape)
        mm[:] = arr
        mm_arrays.append(mm)
    return mm


if __name__ == '__main__':
    import sys
    wt_folder = sys.argv[1]
    mut_folder = sys.argv[2]
    out_dir = sys.argv[3]
    run(wt_folder, mut_folder, out_dir)