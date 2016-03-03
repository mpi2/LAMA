#!/usr/bin/env python

import SimpleITK as sitk
import numpy as np
import sys


im_path = sys.argv[1]
im = sitk.ReadImage(im_path)
arr = sitk.GetArrayFromImage(im)

print "\n"
print "sitk Image: im"
print "numpy array: arr"
print "\n"
print "size: {} ".format(arr.size)
print "dims (z,y,x): {}".format(arr.shape)
print "min: {}".format(arr.min())
print "max: {}".format(arr.max())
print "mean: {}".format(np.mean(arr))
print "spacing: {}".format(im.GetSpacing())
print "origin: {}".format(im.GetOrigin())
print "dtype: {}".format(arr.dtype)
