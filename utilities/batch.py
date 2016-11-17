import SimpleITK as sitk
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
import common
from os.path import join, basename, splitext


indir = '/home/neil/work_ssd/tfce_permutations/inputs_3x_shrunk/all'
outdir = '/home/neil/work_ssd/tfce_permutations/inputs_3x_shrunk/all_nii'

paths = common.GetFilePaths(indir)

for path in paths:
    img = sitk.ReadImage(path)
    bn = splitext(basename(path))[0]
    out_path = join(outdir, bn + '.nii')
    sitk.WriteImage(img, out_path)