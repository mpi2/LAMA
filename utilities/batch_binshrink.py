#!/usr/bin/env python

import SimpleITK as sitk
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import common

indir = sys.argv[1]
outdir = sys.argv[2]
scale_factor = int(sys.argv[3])

images = common.GetFilePaths(indir)

for img in images:

    vol = sitk.ReadImage(img)
    shrunk = sitk.BinShrink(vol, (scale_factor, scale_factor, scale_factor))
    bn = os.path.basename(img)
    outpath = os.path.join(outdir, bn)
    sitk.WriteImage(shrunk, outpath)

