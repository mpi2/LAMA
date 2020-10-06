import SimpleITK as sitk
import os, sys
from pathlib import Path
import numpy as np
import scipy.ndimage as ndimage

# All methods for Harwell are in the PEP8 python format. All ITK methods are in CamelCase.
# Contains methods used for masking
# Also, contains general segmentation methods

# #####################################################################################################################
# Masking methods
######################################################################################################################

def auto_thres_mask(option,
                    img,
                    output,
                    bins=False,
                    tight=False,
                    external=False,
                    external2=False,
                    dilate = True):
    """
    Performs an auto-threshold mask

    :param str option: otsu, huang, li, triangle or isodata
    :param obj img: SimpleITK original image object
    :param str output: Output folder
    :param int bins: The number of bins used if a different number to default simpleITK. [Default False]
    :param boolean tight: True if only a small amount of dilation is required. [Default False]
    :param boolean external: If external hole filling (2D hole filling) is required. [Default False]
    :return obj seg: SimpleITK image object of cleaned up wateshed masked binary object
    """
    print("================== Starting auto-threshold ==================")
    print("Auto method:", option)
    print("doing blur")
    seg = sitk.DiscreteGaussian(img, 3)

    print("doing threshold")

    if not bins:
        bins = 128
    seg = sitk.OtsuThreshold(seg, 0, 255, bins)


    seg = clean_up(img,
                   seg,
                   output,
                   option,
                   force_background_clean = False,
                   dilate = dilate,
                   tight = tight,
                   external = external,
                   external2=external2)


    return seg

def clean_up(img,
             seg,
             output,
             name,
             dilate=True,
             force_background_clean=False,
             tight=False,
             external=False,
             external2=False):
    """ Performs a series of standard operations which are required after a masking operation

    These include

        1. Get largest component in the image
        2. Dilate the component
        3. Fill holes
        4. Create a masked original image
        5. Create an overlay of the binary mask onto the original image

    :param obj img: Original simpleITK image object
    :param obj seg: Binary threshold image
    :param str output: Output directory
    :param str name: Name to save to output file names
    :param boolean dilate: If True dilation occurs if false does not [Default True]
    :param boolean force_background_clean: For watershed. Sometimes the main image is the not first component.
                                            [Default False]
    :param boolean tight: True if image to have small amount of dilation [Default False]
    :param boolean external: True if external holes to be filled [Default False]
    :return simpleitk image object;
    """
    print("================== Clean up ==================")
    print("-Extracting main component")
    # Sometimes the background doesnt get assigned to 0
    if force_background_clean:
        # To be sure use 3, get more accurate result if use a lower numbe though.
        seg = seg > 3

    seg = sitk.ConnectedComponent(seg != seg[0, 0, 0])
    #sitk.WriteImage(seg, os.path.join(output, name + "_all_connected.nrrd"))
    seg = sitk.RelabelComponent(seg)
    seg = seg == 1
    #sitk.WriteImage(seg, os.path.join(output, name + "_largest_connected.nrrd"))
    rad1 = int((seg.GetWidth() * 0.05))

    if dilate:
        print("-Dilating")
        print(tight)

        kernal = [rad1, rad1, rad1]
        if tight:
            print("kernal rad: 1")
            kernal = [1, 1, 1]
        else:
            print("kernal rad", rad1)


        dilate = sitk.BinaryDilateImageFilter()
        dilate.SetKernelRadius(kernal)
        dilate.SetKernelType(sitk.sitkBall)
        seg = dilate.Execute(seg)

    #sitk.WriteImage(seg, os.path.join(output, name + "_dilated.nrrd"))

    if external:
        print("-External hole filling")
        seg = fill_image_all(seg)  # NH:  Why is the same call repeated twice?
        seg = fill_image_all(seg)
    else:
        print("-Sitk fill hole")
        seg = sitk.BinaryFillhole(seg)

    # This is where it fails
    # sitk.WriteImage(seg, os.path.join(output, "wholemask_" + sitk.ImageFileReader.GetFileName(img)))
    output = Path(output)
    outpath = output / f"{output.name}.nrrd"
    sitk.WriteImage(seg, str(outpath))
    #print("-Masking")
    #mask = sitk.Mask(img, seg)
    #sitk.WriteImage(mask, os.path.join(output, name + "_masked.nrrd"))

    #print("-Overlay")
    #tif = os.path.join(output, name + "_overlay.tif")
    #if os.path.exists(tif):
    #    os.remove(tif)
    #overlay = sitk.LabelOverlay(sitk.Cast(sitk.RescaleIntensity(img), sitk.sitkUInt8), seg)
    #sitk.WriteImage(overlay, os.path.join(output, name + "_overlay.tif"))

    print("-Calculating voxels")
    stats = sitk.StatisticsImageFilter()
    stats.Execute(seg)
    print("Total voxels of seg:")
    #should we just return stats.GetSum into a csv file
    
    print(stats.GetSum())

    print("finished")

    return seg

def fill_image_all(input):
    """ Fill an image using 2D hole filling (external hole filling)

    Does this in 3 different orientations

    :param obj input: Simpleitk image object
    :return obj filled: SimpleITK image
    """
    # Operation is peformed using scipy so needs to be a numpy array
    npa = sitk.GetArrayFromImage(input)

    print("fill holes in first orientation")
    npa_hole_filled = fill_image(npa)

    print("fill holes in second orientation")
    npa_hole_filled = fill_image(npa_hole_filled, roll=1)

    print("fill holes in third orientation")
    npa_hole_filled = fill_image(npa_hole_filled, roll=0)

    # Need to turn the image back into its original orientation
    transposed = np.transpose(npa_hole_filled, axes=(0, 2, 1))

    # Turn np array to image
    filled = sitk.GetImageFromArray(transposed)

    # make into a savable format
    filled = sitk.Cast(filled, sitk.sitkUInt8)
    return filled


def fill_image(npa, roll=0):
    """ Binary hole fill in any orientation

    Go slice by slice and fill any holes of a 3D image. Orientation is chosen by the roll parameter.

    :param array npa: Numpy array of the image to be processed
    :param int roll: Which axis to image is to be "rolled" in. If use none, 0 and 1. Should do all axis for a 3D image
                    [Default 0]
    :return array npa_hole_filled: Numpy array of image
    """

    npa_hole_filled = None

    # Change orientation
    if roll:
        loop = np.rollaxis(npa, roll)
    else:
        loop = npa

    # loop through each slice
    for slice_ in loop:
        slice_fill = ndimage.binary_fill_holes(slice_).astype(int)
        if npa_hole_filled is None:
            npa_hole_filled = slice_fill
        else:
            npa_hole_filled = np.dstack((npa_hole_filled, slice_fill))

    return npa_hole_filled
