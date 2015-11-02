import SimpleITK as sitk
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import harwellimglib as hil

indir = sys.argv[1]
outpath = sys.argv[2]


images = hil.GetFilePaths(indir)

#sum all images together
first_img = sitk.ReadImage(images[0])
summed = sitk.GetArrayFromImage(first_img)

summed = summed.astype(np.float32)
for image in images[1:]:#Ommit the first as we have that already
    itk_img = sitk.ReadImage(image)
    np_array = sitk.GetArrayFromImage(itk_img)
    try:
        summed = summed + np_array
    except ValueError as e:
        print "Numpy can't average this volume {0}".format(image)

#Now make average. Do it in numpy as I know how
avg_array = summed/len(images)
avg_img = sitk.GetImageFromArray(avg_array)
print avg_array.max(), avg_array.min(), avg_array.mean()

sitk.WriteImage(avg_img, outpath)

