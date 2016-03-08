#!/usr/bin/env python

import SimpleITK as sitk
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import common

indir = sys.argv[1]
outdir = sys.argv[2]
origin = float(sys.argv[3])
spacing = float(sys.argv[4])

images = common.GetFilePaths(indir)

for img in images:

    vol = sitk.ReadImage(os.path.join(indir, img))
    vol.SetOrigin((origin, origin, origin))
    vol.SetSpacing((spacing, spacing, spacing))
    bn = os.path.basename(img)
    outpath = os.path.join(outdir, bn)
    sitk.WriteImage(vol, outpath, True)

