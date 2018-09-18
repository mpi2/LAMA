#!/usr/bin/env python

import os
import sys
from os.path import join
import numpy as np
import common
from collections import OrderedDict
import tempfile
from logzero import logger as logging
from itertools import chain

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
    given paths to registered images, apply linear normalisation so that the mean of the roi across all images are
    the same.

    Create new diretories and place the normalised images in

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
        numpy.ndarray: use a mask where values == 1 to generate value from

    Returns
    -------
    tuple
        paths to wt normalised images (wts, muts)
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

    all_roi_values = []

    shape = list(memmap_wt_imgs.values())[0].shape

    if isinstance(roi, np.ndarray):
        logging.info('Normalising images to mask')
        roi_mask = roi == 1

    elif isinstance(roi, common.Roi):
        logging.info('Normalising images to roi:  x1:{}, y1:{}, z1:{}  x2:{}, '
                     'y2:{}, z2:{}'.format(roi.x1, roi.x2, roi.y1, roi.y2, roi.z1, roi.z2))
        roi_mask = np.s_[roi.z1: roi.z2, roi.y1: roi.y2, roi.x1: roi.x2]

    # Get mean value across all images
    for basename, imgarr in dict(memmap_wt_imgs, **memmap_mut_imgs).items():  # ** unpacks the second dict
        if isinstance(roi_mask, np.ndarray):  # mask is 1D?
            imgarr = imgarr.ravel()
        all_roi_values.extend(list(imgarr[roi_mask]))
    mean_intensity_all_specimens = np.mean(all_roi_values)

    wt_norm_paths = _apply_linear_normalisation(memmap_wt_imgs, mean_intensity_all_specimens, wt_out_dir, roi_mask, shape, 'wt')
    mut_norm_paths = _apply_linear_normalisation(memmap_mut_imgs, mean_intensity_all_specimens, mut_out_dir, roi_mask, shape, 'mut')

    return wt_norm_paths, mut_norm_paths


def _apply_linear_normalisation(arrays, mean, out_dir, roi_mask, shape, name):

    normalised_paths = []
    log_bits = ["#### Normalisation values for {}".format(name)]

    for basename, imgarr in arrays.items():
        if isinstance(roi_mask, np.ndarray):
            imgarr = imgarr.ravel()
        try:
            roi_values = imgarr[roi_mask]
        except IndexError as e:
            raise

        specimen_mean_roi = roi_values.mean()
        # find deviation from average
        meandiff = specimen_mean_roi - mean
        log_bits.append("{}:    {}".format(basename, str(meandiff)))

        # subtract the difference from each pixel
        try:
            imgarr -= meandiff.astype(np.uint16)  # imagarr = 16bit meandiff = 64bit
        except TypeError:  # Could be caused by imgarr being a short
            imgarr -= int(np.round(meandiff))

        outpath = os.path.join(out_dir, basename)
        common.write_array(imgarr.reshape(shape), outpath)
        normalised_paths.append(outpath)

    log_bits.append('####')
    logging.info("\n".join(log_bits))

    return normalised_paths


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
