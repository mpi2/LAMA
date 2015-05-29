#!/usr/bin/env python

import SimpleITK as sitk
import sys
import os
import numpy as np
import harwellimglib as hil



def normalise(indir, outdir, z, y, x):
    print 'normming', indir
    outdir = os.path.abspath(outdir)
    indir = os.path.abspath(indir)

    for path in hil.GetFilePaths(indir):
        try:
            im = sitk.ReadImage(path)
        except Exception as e:
            print "cant read {}: {}".format(path, e)

        arr = sitk.GetArrayFromImage(im)
        roi = arr[z[0]: z[1], y[0]: y[1], x[0]: x[1]]

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
