#!/usr/bin/env python

import os
import sys
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


def memorymap_data(file_paths):
    imgs = OrderedDict()
    for imgpath in file_paths:
        basename = os.path.basename(imgpath)
        loader = common.LoadImage(imgpath)
        if not loader:
            logging.error("Problem normalising image: {}".format(loader.error_msg))
            sys.exit()
        arr = loader.array
        t = tempfile.TemporaryFile()
        m = np.memmap(t, dtype=arr.dtype, mode='w+', shape=arr.shape)
        m[:] = arr
        imgs[basename] = m
    return imgs


def normalise(wt_paths, mut_paths, outdir, start_indices, end_indices):
    """
    Given a set of images normalise to the mean of the roi and save to disk
    :param indir:
    :param outdir:
    :param start_indexes: iterable, start indices of roi [x, y, z]
    :param lengths: iterable, lengths of roi dimension [x, y, z]
    :return:
    """

    if not isinstance(wt_paths, list) and os.path.isdir(wt_paths):
        wt_paths = common.GetFilePaths(os.path.abspath(wt_paths))
    if not isinstance(mut_paths, list) and os.path.isdir(mut_paths):
        mut_paths = common.GetFilePaths(os.path.abspath(mut_paths))
    outdir = os.path.abspath(outdir)
    common.mkdir_force(outdir)
    wt_out_dir = join(outdir, "wild_type")
    common.mkdir_force(wt_out_dir)
    mut_out_dir = join(outdir, "mutant")
    common.mkdir_force(mut_out_dir)

    memmap_wt_imgs = memorymap_data(wt_paths)
    memmap_mut_imgs = memorymap_data(mut_paths)
    # xyz, from config file
    xs, ys, zs, = start_indices
    xe, ye, ze, = end_indices

    wt_norms_paths = []
    mut_norm_paths = []
    wt_means = []
    mut_means = []
    for basename, imgarr in memmap_wt_imgs.iteritems():
        wt_means.extend(list(imgarr[zs: ze, ys: ye, xs: xe]))
    for basename, imgarr in memmap_mut_imgs.iteritems():
        mut_means.extend(list(imgarr[zs: ze, ys: ye, xs: xe]))
    mean_roi_all = np.mean(wt_means + mut_means)

    logging.info('Normalising images to roi:  x1:{}, y1:{}, z1:{}  x2:{}, '
                 'y2:{}, z2:{} to a mean of {}'.format(xs, ys, zs, xe, ye, ze, mean_roi_all))

    for basename, imgarr in memmap_wt_imgs.iteritems():
        roi = imgarr[zs: ze, ys: ye, xs: xe]   # region of interest, as outlined by command line args

        # Normalise
        # mean of the region of interest
        meanroi = roi.mean()
        # find deviation from average
        meandiff = meanroi - mean_roi_all
        # subtract the difference from each pixel
        try:
            imgarr -= meandiff
        except TypeError:  # Could be caused by imgarr being a short
            imgarr -= int(np.round(meandiff))

        outpath = os.path.join(wt_out_dir, basename)
        common.write_array(imgarr, outpath)
        wt_norms_paths.append(outpath)

    for basename, imgarr in memmap_mut_imgs.iteritems():
        roi = imgarr[zs: ze, ys: ye, xs: xe]  # region of interest, as outlined by command line args

        # Normalise
        meanroi = roi.mean()  # mean of the region of interest
        meandiff = meanroi - mean_roi_all  # finds deviation from reference
        try:
            imgarr -= meandiff
        except TypeError:  # Could be caused by imgarr being a short
            imgarr -= int(np.round(meandiff))

        outpath = os.path.join(mut_out_dir, basename)
        common.write_array(imgarr, outpath)
        mut_norm_paths.append(outpath)

    return wt_norms_paths, mut_norm_paths

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-w', dest='wt', help='wt dir', required=True)
    parser.add_argument('-m', dest='mut', help='mut dir', required=True)
    parser.add_argument('-o', dest='output', help='output dir', required=True)
    parser.add_argument('-s', dest='starts', help='start indices (x, y, z)', required=True, nargs=3, type=int)
    parser.add_argument('-e', dest='ends', help='end indices (x, y, z)', required=True, nargs=3, type=int)

    args = parser.parse_args()
    wt = common.GetFilePaths(args.wt)

    mut = common.GetFilePaths(args.mut)
    normalise(wt, mut, args.output, args.starts, args.ends)
