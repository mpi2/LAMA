#!/usr/bin/env python

import SimpleITK as sitk
import os
import numpy as np
import common
from collections import OrderedDict
import tempfile


def memorymap_data(dirs):
    imgs = OrderedDict()
    for d in dirs:
        for imgpath in common.GetFilePaths(d):
            basename = os.path.basename(imgpath)
            arr = common.img_path_to_array(imgpath)
            t = tempfile.TemporaryFile()
            m = np.memmap(t, dtype=arr.dtype, mode='w+', shape=arr.shape)
            m[:] = arr
            imgs[basename] = m
    return imgs


def normalise(indir, outdir, start_indices, end_indices):
    """
    :param indir:
    :param outdir:
    :param start_indexes: iterable, start indices of roi [z, y, x]
    :param lengths: iterable, lengths of roi dimension [z, y, x]
    :return:
    """

    outdir = os.path.abspath(outdir)
    indir = os.path.abspath(indir)

    memmap_imgs = memorymap_data([indir])
    zs, ys, xs = start_indices
    ze, ye, xe = end_indices

    means = []
    for basename, imgarr in memmap_imgs.iteritems():
        means.extend(list(imgarr[zs: ze, ys: ye, xs: xe]))
    mean_roi_all = np.mean(means)

    for basename, imgarr in memmap_imgs.iteritems():
        roi = imgarr[zs: ze, ys: ye, xs: xe]   # region of interest, as outlined by command line args

        # Normalise
        meanroi = roi.mean()                # mean of the region of interest
        meandiff = meanroi - mean_roi_all        # finds deviation from reference
        print
        print "meanroi=", meanroi
        print "meanref=", mean_roi_all
        print "meandiff=", meandiff

        imgarr -= meandiff                  # subtracts the difference from each pixel

        outimg = sitk.GetImageFromArray(imgarr)
        outpath = os.path.join(outdir, basename)
        sitk.WriteImage(outimg, outpath)
        print 'done'

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input', help='input dir', required=True)
    parser.add_argument('-o', dest='output', help='output dir', required=True)
    parser.add_argument('-s', dest='starts', help='start indices (z, y, x)', required=True, nargs=3, type=int)
    parser.add_argument('-e', dest='ends', help='end indices (z, y, x)', required=True, nargs=3, type=int)

    args = parser.parse_args()

    normalise(args.input, args.output, args.starts, args.ends)
