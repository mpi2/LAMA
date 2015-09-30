"""
Given a subarray cube, and the orginal dimensions of that subarray, insert it back into a vactor image of the same
dimensiuons of the one it cam from
"""

import numpy as np
import SimpleITK as sitk

gauss_width = 2
gaus_var = 0.5

shrunk_array_dims = ((177, 205), (69, 91), (79, 99))

image_to_warp_path = '/home/neil/work/defined_abnormalites/test_data/20150411_RAPSN_E14.5_16.1f_WT_XY_rec_scaled_4.6823_pixel_14.0001.nrrd'
sub_array_path = '/home/neil/work/defined_abnormalites/from_henrik/def_1.nrrd'
warp_out = '/home/neil/work/defined_abnormalites/from_henrik/warp-out/warped.nrrd'

original_image = sitk.ReadImage(image_to_warp_path)
original_array = sitk.GetArrayFromImage(original_image)
original_dims = original_array.shape

shrunk_array = sitk.GetArrayFromImage(sitk.ReadImage(sub_array_path))

zero_def_field = np.zeros(list(original_dims) + [3])

# Stik the shrunk organ eformation field into the zero vector array

s = shrunk_array_dims

zero_def_field[s[0][0]:s[0][1], s[1][0]: s[1][1], s[2][0]: s[2][1]] = shrunk_array

# Invert the deformation field or we get exapnsion with Warp()

zero_def_field = zero_def_field * -1

# Blur the vectors to get a smoother result

zvectors = sitk.GetImageFromArray(zero_def_field[:, :, :, 0])
yvectors = sitk.GetImageFromArray(zero_def_field[:, :, :, 1])
xvectors = sitk.GetImageFromArray(zero_def_field[:, :, :, 2])


zblured = sitk.DiscreteGaussian(zvectors, variance=1.0, maximumKernelWidth=gauss_width, maximumError=0.01, useImageSpacing=False)
yblured = sitk.DiscreteGaussian(yvectors, variance=1.0, maximumKernelWidth=gauss_width, maximumError=0.01, useImageSpacing=False)
xblured = sitk.DiscreteGaussian(xvectors, variance=1.0, maximumKernelWidth=gauss_width, maximumError=0.01, useImageSpacing=False)


# reintroduce the blurred vectors back into the vector field

zero_def_field[:, :, :, 0] = sitk.GetArrayFromImage(zblured)
zero_def_field[:, :, :, 1] = sitk.GetArrayFromImage(yblured)
zero_def_field[:, :, :, 2] = sitk.GetArrayFromImage(xblured)

# Warp the original image
deformation_img = sitk.GetImageFromArray(zero_def_field)

warped = sitk.Warp(original_image, deformation_img)

sitk.WriteImage(warped, warp_out)





