#! /usr/bin/env python

import os
import numpy as np
import SimpleITK as sitk


def run(vol_paths, outdir):

    for path in vol_paths:
        try:
            img = sitk.ReadImage(path)
        except RuntimeError:
            # Not an image
            continue

        arr = sitk.GetArrayFromImage(img).astype(np.float16)
        vector_magnitudes = np.sqrt((arr*arr).sum(axis=3))
        mag_img = sitk.GetImageFromArray(vector_magnitudes)
        base = os.path.basename(path)
        outpath = os.path.join(outdir, base)
        sitk.WriteImage(mag_img, outpath, True)


if __name__ == '__main__':
    import sys


    indir = sys.argv[1]
    outdir = sys.argv[2]
    vols = []
    for path in os.listdir(indir):
        vols.append(os.path.join(indir, path))
    run(vols, outdir)
