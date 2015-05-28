#!/usr/bin/env python

import SimpleITK as sitk
import sys
import os
import numpy as np
import harwellimglib as hil



def normalise(indir, outdir, size, index):
    print 'normming'
    outdir = os.path.abspath(outdir)
    indir = os.path.abspath(indir)

    for path in hil.GetFilePaths(indir):
        try:
            im = sitk.ReadImage(path)
        except:
            print "cant read {}".format(path)

        roi = sitk.RegionOfInterest(im, size, index)
        arr = sitk.GetArrayFromImage(roi)

        mean_ = float(np.mean(arr))

        print 'mean', mean_, type(mean_)

        #cast to 8bit
        cast = sitk.Cast(sitk.RescaleIntensity(im), sitk.sitkUInt8)
        # Normalise
        norm = sitk.IntensityWindowing(im, mean_)

        basename = os.path.basename(path)
        outpath = os.path.join(outdir, basename)
        sitk.WriteImage(norm, outpath)


if __name__ == '__main__':
    in_ = sys.argv[1]
    out = sys.argv[2]

    s = [9, 7, 8]
    i = [103, 131, 323]

    normalise(in_, out, s, i)
