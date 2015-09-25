#!/usr/bin/env python

import sys
import SimpleITK as sitk
import numpy as np
import subprocess as sub


deformation_p = sys.argv[1]
img_p = sys.argv[2]
out_path = sys.argv[3]


img = sitk.ReadImage(img_p)
arr = sitk.GetArrayFromImage(img)
shape_str = [str(x) for x in arr.shape]
shape =  arr.shape

deformation = np.fromfile(deformation_p, dtype=np.float32)

try:
    sub.check_output(['WarpImage3D', '-I', img_p, '-d', deformation_p, '-r', shape_str[1], '-c', shape_str[2], '-s', shape_str[0],
                  '-o', out_path])
except Exception:
    pass

deformed = np.fromfile(out_path, dtype=np.byte).reshape(shape)
im_out = sitk.GetImageFromArray(deformed)
sitk.WriteImage(im_out, 'final.nrrd')
