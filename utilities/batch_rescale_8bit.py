batch_rescale_8bit.py#!/usr/bin/env python

import SimpleITK as sitk
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
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
    cast = sitk.Cast(sitk.RescaleIntensity(img), sitk.sitkUInt8)
    basename = os.path.basename(path)
    outpath = os.path.join(outdir, basename)
    sitk.WriteImage(cast, outpath)


