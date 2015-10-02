"""
Given a subarray cube, and the orginal dimensions of that subarray relative to the original image, insert it into a
vector image of the same dimensiuons of the one it cam from
"""

import numpy as np
import SimpleITK as sitk
import scipy.ndimage as ndimage


gauss_width = 4
gaus_var = 1.0
def_mask_dilate_val = 4  # Dilation value for the vector field mask

shrunk_array_dims = ((125, 151), (136, 167), (122, 153))

image_to_warp_path = '/home/neil/work/defined_abnormalites/registration/synthetic_mutants/20141127_PRKAB2_E14.5_15.1G_WT_XX/20141127_PRKAB2_E14.5_15.1G_WT_XX_rec_scaled_4.6878_pixel_14.tif'
def_roi_path = '/home/neil/work/defined_abnormalites/registration/synthetic_mutants/20141127_PRKAB2_E14.5_15.1G_WT_XX/shrunk/def_1.nrrd'
warp_out = '/home/neil/work/defined_abnormalites/registration/synthetic_mutants/20141127_PRKAB2_E14.5_15.1G_WT_XX/warped/warped.nrrd'
def_field_out_path = '/home/neil/work/defined_abnormalites/registration/synthetic_mutants/20141127_PRKAB2_E14.5_15.1G_WT_XX/warped/roi_def_field.mhd'

# This is the original mask. Dilate it and use it as a vector mask
roi_mask_path = '/home/neil/work/defined_abnormalites/registration/synthetic_mutants/20141127_PRKAB2_E14.5_15.1G_WT_XX/20141127_PRKAB2_E14.5_15.1G_WT_XX_adrenal_roi.nrrd'


def blur(array):
    """
    Try and blut the directions
    """
    # Blur the vectors to get a smoother result

    zvectors = sitk.GetImageFromArray(array[:, :, :, 0])
    yvectors = sitk.GetImageFromArray(array[:, :, :, 1])
    xvectors = sitk.GetImageFromArray(array[:, :, :, 2])

    zblured = sitk.DiscreteGaussian(zvectors, variance=gaus_var, maximumKernelWidth=gauss_width, maximumError=0.01, useImageSpacing=False)
    yblured = sitk.DiscreteGaussian(yvectors, variance=gaus_var, maximumKernelWidth=gauss_width, maximumError=0.01, useImageSpacing=False)
    xblured = sitk.DiscreteGaussian(xvectors, variance=gaus_var, maximumKernelWidth=gauss_width, maximumError=0.01, useImageSpacing=False)

    array[:, :, :, 0] = sitk.GetArrayFromImage(zblured)
    array[:, :, :, 1] = sitk.GetArrayFromImage(yblured)
    array[:, :, :, 2] = sitk.GetArrayFromImage(xblured)

    return array


def mask_deformation_field(mask_path, deformation_field):
    mask_img = sitk.ReadImage(mask_path)
    dilated_mask_img = sitk.BinaryDilate(mask_img, def_mask_dilate_val)
    dilated_mask_array = sitk.GetArrayFromImage(dilated_mask_img)

    deformation_field[dilated_mask_array == 0] = 0.0
    return deformation_field


original_image = sitk.ReadImage(image_to_warp_path)
original_array = sitk.GetArrayFromImage(original_image)
original_dims = original_array.shape

shrunk_array = sitk.GetArrayFromImage(sitk.ReadImage(def_roi_path))
masked_shrunk_array = mask_deformation_field(roi_mask_path, shrunk_array)

zero_vector_field = np.zeros(list(original_dims) + [3])

# Stick the shrunk organ deformation field into the zero vector array
s = shrunk_array_dims
zero_vector_field[s[0][0]:s[0][1], s[1][0]: s[1][1], s[2][0]: s[2][1]] = masked_shrunk_array

blurred_array = blur(zero_vector_field)

# Invert the deformation field or we get exapnsion with Warp()
blurred_array = blurred_array * -1

# Warp the original image
deformation_img = sitk.GetImageFromArray(blurred_array)
warped = sitk.Warp(original_image, deformation_img, interpolator=sitk.sitkBSpline)

sitk.WriteImage(warped, warp_out)

# Write to file the ROI def field for visual analysis
def_field_roi = blurred_array[s[0][0]:s[0][1], s[1][0]: s[1][1], s[2][0]: s[2][1]]
sitk.WriteImage(sitk.GetImageFromArray(def_field_roi), def_field_out_path)





