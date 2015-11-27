#!/usr/bin/env python

"""
Given a subarray cube, and the orginal dimensions of that subarray relative to the original image, insert it into a
vector image of the same dimensiuons of the one it cam from
"""

import numpy as np
import SimpleITK as sitk
import scipy.ndimage as ndimage


gauss_width = 4
gaus_var = 1.0
def_mask_dilate_val = 2  # Dilation value for the vector field mask


# This is the original mask. Dilate it and use it as a vector mask
roi_mask_path = '/home/neil/work/defined_abnormalites/registration/synthetic_mutants/20141127_PRKAB2_E14.5_15.1G_WT_XX/20141127_PRKAB2_E14.5_15.1G_WT_XX_adrenal_roi.nrrd'


def blur(def_img):
    """
    Try and blut the directions
    """
    # Blur the vectors to get a smoother result

    array = sitk.GetArrayFromImage(def_img)

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


def mask_deformation_field(mask_path, deformation_field, outpath):
    mask_img = sitk.ReadImage(mask_path)
    dilated_mask_img = sitk.BinaryDilate(mask_img, def_mask_dilate_val)
    dilated_mask_array = sitk.GetArrayFromImage(dilated_mask_img)

    deformation_field[dilated_mask_array == 0] = 0.0
    return deformation_field


def run(image_to_warp_path, def_field_path, outpath):

    vol_to_warp = sitk.ReadImage(image_to_warp_path)
    def_field = sitk.ReadImage(def_field_path)
    blurred_def_array = blur(def_field)

    # Invert the deformation field or we get exapnsion with Warp()
    blurred_def_array= blurred_def_array * -1

    # Warp the original image
    deformation_img = sitk.GetImageFromArray(blurred_def_array)
    warped = sitk.Warp(vol_to_warp, deformation_img, interpolator=sitk.sitkBSpline)

    sitk.WriteImage(warped, outpath)




import argparse


parser = argparse.ArgumentParser("Apply deformation field to image")
parser.add_argument('-v', '--vol', dest='vol', help='volume to applpy def field to', required=True)
parser.add_argument('-d', '--def', dest='def_field', help='deformation field', required=True)
parser.add_argument('-o', '--out', dest='out', help='out file path with extension', required=True)

#parser.add_argument('-n', '--label_names', dest='label_names', help='csv files with label names', required=True)
args = parser.parse_args()
run(args.vol, args.def_field, args.out)
