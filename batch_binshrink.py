#!/usr/bin/env python

import SimpleITK as sitk
import sys
import os
sys.path.append("..")
import harwellimglib as hil

FACTOR = 2

indir = sys.argv[1]
outdir = sys.argv[2]

images = hil.GetFilePaths(indir)

for img in images:
    vol = sitk.ReadImage(img)
    shrunk = sitk.BinShrink(vol, (FACTOR, FACTOR, FACTOR))
    bn = os.path.basename(img)
    outpath = os.path.join(outdir, bn)
    sitk.WriteImage(shrunk, outpath)



