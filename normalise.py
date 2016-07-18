#!/usr/bin/env python

import os
from os.path import join
import numpy as np
import common
from collections import OrderedDict
import tempfile
import logging

try:
    from skimage.draw import line_aa
    skimage_available = True
except ImportError:
    skimage_available = False


def memorymap_data(dirs):
    imgs = OrderedDict()
    for d in dirs:
        for imgpath in common.GetFilePaths(d):
            basename = os.path.basename(imgpath)
            loader = common.LoadImage(imgpath)
            if not loader:
                logging.error("Problem normalising image: {}".format(loader.error_msg))
            arr = loader.array
            t = tempfile.TemporaryFile()
            m = np.memmap(t, dtype=arr.dtype, mode='w+', shape=arr.shape)
            m[:] = arr
            imgs[basename] = m
    return imgs


def normalise(wt_dir, mut_dir, outdir, start_indices, end_indices):
    """
    Given a set of images normalise to the mean of the roi and save to disk
    :param indir:
    :param outdir:
    :param start_indexes: iterable, start indices of roi [x, y, z]
    :param lengths: iterable, lengths of roi dimension [x, y, z]
    :return:
    """

    outdir = os.path.abspath(outdir)
    common.mkdir_force(outdir)
    wt_out_dir = join(outdir, "willd_type")
    common.mkdir_force(wt_out_dir)
    mut_out_dir = join(outdir, "mutant")
    common.mkdir_force(mut_out_dir)
    wt_indir = os.path.abspath(wt_dir)
    mut_indir =os.path.abspath(mut_dir)

    memmap_wt_imgs = memorymap_data([wt_indir])
    memmap_mut_imgs = memorymap_data([mut_indir])
    # xyz, from config file
    xs, ys, zs, = start_indices
    xe, ye, ze, = end_indices

    means = []

    wt_norms_paths = []
    mut_norm_paths = []
    for basename, imgarr in memmap_wt_imgs.iteritems():
        means.extend(list(imgarr[zs: ze, ys: ye, xs: xe]))
    for basename, imgarr in memmap_mut_imgs.iteritems():
        means.extend(list(imgarr[zs: ze, ys: ye, xs: xe]))
    mean_roi_all = np.mean(means)

    for basename, imgarr in memmap_wt_imgs.iteritems():
        roi = imgarr[zs: ze, ys: ye, xs: xe]   # region of interest, as outlined by command line args

        # Normalise
        meanroi = roi.mean()                # mean of the region of interest
        meandiff = meanroi - mean_roi_all        # finds deviation from reference
        imgarr -= meandiff                  # subtracts the difference from each pixel

        outpath = os.path.join(wt_out_dir, basename)
        common.write_array(imgarr, outpath)
        wt_norms_paths.append(outpath)
    print 'normalising wild types done'

    for basename, imgarr in memmap_mut_imgs.iteritems():
        roi = imgarr[zs: ze, ys: ye, xs: xe]  # region of interest, as outlined by command line args

        # Normalise
        meanroi = roi.mean()  # mean of the region of interest
        meandiff = meanroi - mean_roi_all  # finds deviation from reference
        imgarr -= meandiff  # subtracts the difference from each pixel

        outpath = os.path.join(mut_out_dir, basename)
        common.write_array(imgarr, outpath)
        mut_norm_paths.append(outpath)
    print 'normalising mutants done'
    return wt_norms_paths, mut_norm_paths

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input', help='input dir', required=True)
    parser.add_argument('-o', dest='output', help='output dir', required=True)
    parser.add_argument('-s', dest='starts', help='start indices (x, y, z)', required=True, nargs=3, type=int)
    parser.add_argument('-e', dest='ends', help='end indices (x, y, z)', required=True, nargs=3, type=int)

    args = parser.parse_args()

    normalise(args.input, args.output, args.starts, args.ends)
