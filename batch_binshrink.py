#!/usr/bin/env python

import SimpleITK as sitk
import sys
import os
sys.path.append("..")
import harwellimglib as hil



indir = sys.argv[1]
outdir = sys.argv[2]
shrink_factor = int(sys.argv[3])

images = hil.GetFilePaths(indir)

for img in images:
    vol = sitk.ReadImage(img)
    shrunk = sitk.BinShrink(vol, (shrink_factor, shrink_factor, shrink_factor))
    bn = os.path.basename(img)
    outpath = os.path.join(outdir, bn)
    sitk.WriteImage(shrunk, outpath)



