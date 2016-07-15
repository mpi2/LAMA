#!/usr/bin/env python

import SimpleITK as sitk
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
import common


if len(sys.argv) > 2:
    indir = sys.argv[1]
    outdir = sys.argv[2]
else:
    indir = sys.argv[1]
    # Relplace volumes
    outdir = sys.argv[1]


paths = common.GetFilePaths(indir)

for path in paths:
    img = sitk.ReadImage(path)
    arr = sitk.GetArrayFromImage(img)
    arr /= 256.0
    arr = arr.astype(np.uint8)
    out_img = sitk.GetImageFromArray(arr)
    basename = os.path.basename(path)
    outpath = os.path.join(outdir, basename)
    sitk.WriteImage(out_img, outpath, True)


