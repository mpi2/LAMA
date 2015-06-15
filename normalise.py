#!/usr/bin/env python

import SimpleITK as sitk
import sys
import os
import numpy as np
import harwellimglib as hil

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

    for path in hil.GetFilePaths(indir):
        try:
            im = sitk.ReadImage(path)
        except Exception as e:
            print "cant read {}: {}".format(path, e)
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
    in_ = sys.argv[1]
    out = sys.argv[2]

    z = (325, 331)
    y = (133, 139)
    x = (104,110)

    normalise(in_, out, z, y, x)
