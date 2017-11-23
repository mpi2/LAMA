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


def normalise(wt_paths, mut_paths, outdir, roi):
    """
    
    Parameters
    ----------
    wt_paths: str
        input paths for normalising
    mut_paths
        input paths for normalising
    outdir: str
        path to put nomalised images
    roi: 
        commom.Roi: roi to generate normalisation value from
        numpy.ndarray: use a mask where values == 1 to gerneate value from
    """

    if not isinstance(wt_paths, list) and os.path.isdir(wt_paths):
        wt_paths = common.get_file_paths(os.path.abspath(wt_paths))
    if not isinstance(mut_paths, list) and os.path.isdir(mut_paths):
        mut_paths = common.get_file_paths(os.path.abspath(mut_paths))
    outdir = os.path.abspath(outdir)
    common.mkdir_force(outdir)
    wt_out_dir = join(outdir, "wild_type")
    common.mkdir_force(wt_out_dir)
    mut_out_dir = join(outdir, "mutant")
    common.mkdir_force(mut_out_dir)

    memmap_wt_imgs = memorymap_data(wt_paths)
    memmap_mut_imgs = memorymap_data(mut_paths)
    # xyz, from config file

    wt_norms_paths = []
    mut_norm_paths = []
    wt_means = []
    mut_means = []

    shape = memmap_wt_imgs.values()[0].shape

    if isinstance(roi, np.ndarray):
        logging.info('Normalising images to mask')
        indx = roi == 1
    elif isinstance(roi, common.Roi):
        logging.info('Normalising images to roi:  x1:{}, y1:{}, z1:{}  x2:{}, '
                     'y2:{}, z2:{}'.format(roi.x1, roi.x2, roi.y1, roi.y2, roi.z1, roi.z2))
        indx = np.s_[roi.z1: roi.z2, roi.y1: roi.y2, roi.x1: roi.x2]

    for basename, imgarr in memmap_wt_imgs.iteritems():
        if isinstance(roi, np.ndarray):
            imgarr = imgarr.ravel()
        wt_means.extend(list(imgarr[indx]))

    for basename, imgarr in memmap_mut_imgs.iteritems():
        if isinstance(roi, np.ndarray):
            imgarr = imgarr.ravel()
        mut_means.extend(list(imgarr[indx]))

    mean_roi_all = np.mean(wt_means + mut_means)

    for basename, imgarr in memmap_wt_imgs.iteritems():
        if isinstance(roi, np.ndarray):
            imgarr = imgarr.ravel()
        try:
            roi_values = imgarr[indx]   # region of interest, as outlined by command line args
        except IndexError as e:
            raise

        # Normalise
        # mean of the region of interest
        meanroi = roi_values.mean()
        # find deviation from average
        meandiff = meanroi - mean_roi_all
        # subtract the difference from each pixel
        try:
            imgarr -= meandiff.astype(np.uint16)  # imagarr = 16bit meandiff = 64bit
        except TypeError:  # Could be caused by imgarr being a short
            imgarr -= int(np.round(meandiff))

        outpath = os.path.join(wt_out_dir, basename)
        common.write_array(imgarr.reshape(shape), outpath)
        wt_norms_paths.append(outpath)

    for basename, imgarr in memmap_mut_imgs.iteritems():
        if isinstance(roi, np.ndarray):
            imgarr = imgarr.ravel()
        try:
            roi_values = imgarr[indx]  # region of interest, as outlined by command line args
        except IndexError:
            raise

        # Normalise
        meanroi = roi_values.mean()  # mean of the region of interest
        meandiff = meanroi - mean_roi_all  # finds deviation from reference
        try:
            imgarr -= meandiff
        except TypeError:  # Could be caused by imgarr being a short
            imgarr -= int(np.round(meandiff))

        outpath = os.path.join(mut_out_dir, basename)
        common.write_array(imgarr.reshape(shape), outpath)
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
    wt = common.get_file_paths(args.wt)

    mut = common.get_file_paths(args.mut)
    normalise(wt, mut, args.output, args.starts, args.ends)
