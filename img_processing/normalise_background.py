#!/usr/bin/env python

import SimpleITK as sitk
import os
import numpy as np
import common


def normalise(indir, outdir, start_indices, end_indices):
    """
    :param indir:
    :param outdir:
    :param start_indexes: iterable, start indices of roi [z, y, x]
    :param lengths: iterable, lengths of roi dimension [z, y, x]
    :return:
    """

    zs, ys, xs = start_indices
    ze, ye, xe = end_indices

    outdir = os.path.abspath(outdir)
    indir = os.path.abspath(indir)

    for path in common.get_file_paths(indir):
        try:
            im = sitk.ReadImage(path)
        except Exception as e:
            print("cant read {}: {}".format(path, e))
            raise

        arr = sitk.GetArrayFromImage(im)
        roi = arr[zs: ze, ys: ye, xs: xe]

        mean_ = float(np.mean(roi))

        # Normalise
        norm = sitk.IntensityWindowing(im, mean_)

        basename = os.path.basename(path)
        outpath = os.path.join(outdir, basename)
        sitk.WriteImage(norm, outpath)


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input', help='input dir', required=True)
    parser.add_argument('-o', dest='output', help='output dir', required=True)
    parser.add_argument('-s', dest='starts', help='start indices (z, y, x)', required=True, nargs=3, type=int)
    parser.add_argument('-e', dest='ends', help='end indices (z, y, x)', required=True, nargs=3, type=int)

    args = parser.parse_args()

    normalise(args.input, args.output, args.starts, args.ends)
