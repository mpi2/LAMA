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
def ws_mask(img,
            output,
            sig=1.25,
            ws_level=0.9,
            force_background=False,
            watershed_line=False,
            fully_connected=False,
            tight=False,
            external=True,
            external2=False,
            dilate = True):

    """
    Perform a watershed mask and clean up

    :param obj img: Original image in SimpleITK image format
    :param str output: Output folder
    :param float sig: Sigma (STDEV) for gradient calculation
    :param float ws_level: Watershed level - Lower the number more components
    :param boolean force_background: Can force the mask to ignore the first few components.
    (Required if another object in background)
    :param boolean watershed_line: If a line is required between each watershed segmentation.
    :param boolean fully_connected: Additional watershed option
    :param boolean tight: If the mask should be tight or not [Default False]
    :param boolean external: If external hole filling (2D hole filling) is required.
    :return obj seg: SimpleITK image object of cleaned up wateshed masked binary object
    """
    print("================== Starting watershed ==================")
    print("Getting gaussian gradient magnitude")
    seg = sitk.GradientMagnitudeRecursiveGaussian(img, sigma=sig)
    sitk.WriteImage(seg, os.path.join(output, "gradient_4_ws.nrrd"))

    print("Doing watershed")
    seg = sitk.MorphologicalWatershed(seg,
                                      level=ws_level,
                                      markWatershedLine=watershed_line,
                                      fullyConnected=fully_connected)

    sitk.WriteImage(seg, os.path.join(output, "ws.nrrd"))

    print("performing clean up")
    seg = clean_up(img,
                   seg,
                   output,
                   "ws",
                   dilate = dilate,
                   force_background_clean = force_background,
                   tight = tight,
                   external = external,
                   external2=external2)

    return seg

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
    if option == "otsu":
        if not bins:
            bins = 128
        seg = sitk.OtsuThreshold(seg, 0, 255, bins)
    elif option == "huang":
        if not bins:
            bins = 128
        seg = sitk.HuangThreshold(seg, 0, 255, bins)
    elif option == "li":
        if not bins:
            bins = 256
        seg = sitk.LiThreshold(seg, 0, 255, bins)
    elif option == "triangle":
        if not bins:
            bins = 256
        seg = sitk.TriangleThreshold(seg, 0, 255, bins)
    elif option == "isodata":
        if not bins:
            bins = 256
        seg = sitk.IsoDataThreshold(seg, 0, 255, bins)

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


def binary_mask(img,
                output,
                lt=False,
                tight=False,
                external=True,
                external2=False,
                dilate = True):

    """  Performs a binary mask where we choose just a lower threshold value

    Also if no-lower threshold use, we calculate a binary threshold using the maximum intensity and variance of the
    the first slice

    Saves files as it goes along.

    :param obj img: SimpleITK image object
    :param str output: Output folder
    :param float/boolean lt: False if to lower-threshold to be calculated automatically, float if manually chosen
                            [Default False]
    :param boolean tight: True if tight fitting mask with only small amount of dilation [Default: False]
    :param boolean external: True if 2D hole filling (external hole filling) required [Default: True]
    """
    print("==================Started binary masking==================")
    # If lower threshold is assigned just perform a normal threshold using this value
    # Simple lower threshold determined binary mask. Upper threshold just the top value
    if lt:
        print("Using manual binary with lower thershold of", lt)
        print("blurring")
        seg = sitk.DiscreteGaussian(img, 1)
        print("thresholding")
        sif = sitk.StatisticsImageFilter()
        sif.Execute(seg)
        max_ = sif.GetMaximum()
        seg = sitk.BinaryThreshold(
            seg,
            lowerThreshold=int(lt),
            upperThreshold=max_,
            insideValue=255,
            outsideValue=0)

        print("saving")
        sitk.WriteImage(seg, os.path.join(output, "binary.nrrd"))
        # Perform clean up stage
        clean_up(img, seg, output, "binary", dilate = dilate, tight = tight, external = external, external2=external2)
        return

    print("Performing binary-auto")
    # Get a region of the first slice to calculate the background intensity

    print("Calculating lower threshold")

    # From the first slice, we get a small subsection which typically doesn't have any other objects. Hence
    # provides a good background slice
    size = list(img.GetSize())
    print(size)
    # x coordinates
    x = int(size[0] * 0.25)
    y = int(size[1] * 0.01)
    print("x,y", x, y)
    size[0] = int(size[0] * 0.25)
    size[1] = int(size[1] * 0.25)
    size[2] = 0
    print("size", size)
    index = [x, y, 0]
    Extractor = sitk.ExtractImageFilter()
    Extractor.SetSize(size)
    Extractor.SetIndex(index)
    slice = Extractor.Execute(img)

    # Get the stats for the extracted region
    stats = sitk.StatisticsImageFilter()
    stats.Execute(slice)
    max = int(stats.GetMaximum())
    mean = int(stats.GetMean())
    sigma = stats.GetSigma()

    # save the region
    # sitk.WriteImage(slice, os.path.join(output, "slice.tif"))

    # Base the mask on max intensity - 2.5*SD (or whatever variance Sigma is referring to)
    factor = 2.5
    lt = stats.GetMaximum() - factor * stats.GetSigma()
    lt = int(lt)

    print("mean", mean)
    print("max", max)
    print("var", sigma)

    # Just in case this lower threshold is too small. We do a check to see if it is the less than the mean * 1.2
    # If so, the maximum intensity of the selected region is used.
    if mean + (mean * 0.2) > lt:
        print("mean still greater than original proposed threshold, just use max")
        lt = int(max)
    else:
        print("Lower threshold", lt)

    # Now that we have a lower threshold we perform the actual image processing.
    print("blurring")
    seg = sitk.DiscreteGaussian(img, 1)
    print("thresholding")
    seg = sitk.BinaryThreshold(seg, lowerThreshold=lt, upperThreshold=255, insideValue=255, outsideValue=0)
    seg = sitk.Cast(seg, sitk.sitkUInt8)
    sitk.WriteImage(seg, os.path.join(output, "binary.nrrd"))
    print("clean up")
    seg = clean_up(img, seg, output, "binary", dilate = dilate, tight= tight, external = external, external2=external2)


def binary_mask_split(img,
                output,
                lt1,
                lt2,
                split,
                tight=False,
                external=False,
                external2=False,
                dilate = True):
    """
    Can split the 3D image in the Z-axis by a user defined value.
    Then use a different lower threshold for each section.

    :param obj img: SimpleITK image object
    :param str output: Output directory
    :param int lt1: Lower treshold for first section of the image
    :param int lt2: Lower treshold for second section of the image
    :param int split: zslice to split the image
    :param boolean tight: True if image to have small amount of dilation [Default False]
    :param boolean external: True if external holes to be filled [Default False]
    """
    #===================================================================
    # First section of image
    #===================================================================
    blur = sitk.DiscreteGaussian(img, 1)
    size = list(img.GetSize())
    z_depth = size[2]
    size[2] = int(split)
    start = [0, 0, 1]
    first = sitk.RegionOfInterest(blur, size, start)

    print("blurring")

    print("thresholding")
    seg_first = sitk.BinaryThreshold(
            first,
            lowerThreshold=int(lt1),
            upperThreshold=255,
            insideValue=255,
            outsideValue=0)

    #===================================================================
    # Second section of the image
    #===================================================================
    z_depth2 = int(z_depth)-int(split)
    size[2] = int(z_depth2)
    start = [0, 0, (int(split)-1)]
    second = sitk.RegionOfInterest(blur, size, start)

    print("blurring")
    seg_second = sitk.DiscreteGaussian(second, 1)
    print("thresholding")
    seg_second = sitk.BinaryThreshold(
            second,
            lowerThreshold=int(lt2),
            upperThreshold=255,
            insideValue=255,
            outsideValue=0)

    # Put back together again. I use numpy as it is not clear how to do this in simple itk
    npa1 = sitk.GetArrayFromImage(seg_first)
    npa2 = sitk.GetArrayFromImage(seg_second)
    npa_all = np.vstack((npa1, npa2))
    # Turn np array to image
    all = sitk.GetImageFromArray(npa_all)

    # make into a savable format
    all = sitk.Cast(all, sitk.sitkUInt8)
    sitk.WriteImage(all , os.path.join(output, "back_together.nrrd"))

    # Perform the clean up
    seg = clean_up(img,
                   all,
                   output,
                   "binary",
                   dilate = dilate,
                   tight = tight,
                   external = external,
                   external2=external2)


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

######################################################################################################################
# Segmentation methods
######################################################################################################################
def zero_crossing(img, output, back=0, fore=255, maximErr=0.1, variance=0.0019):
    '''
    Harwell implementation of ITK ZeroCrossingBasedEdgeDetectionImageFilter
    @param img: sitk image (should be in Float format this is converted in the method)
    @param output: output folder
    @param ZeroCrossingBasedEdgeDetectionImageFilte arguments
    @return out = the edge detected image
    @return edge_ob = the ls object
    '''
    print("***Harwell implementation of ITK ZeroCrossingBasedEdgeDetectionImageFilter***")
    print("Image needs to be a float (now casting)")
    img_f = sitk.Cast(img, sitk.sitkFloat32)
    print("Saving the float image img_f.nrrd")
    sitk.WriteImage(img_f, os.path.join(output, "img_f.nrrd"))
    print("Setting up object")
    edge_ob = sitk.ZeroCrossingBasedEdgeDetectionImageFilter()
    edge_ob.SetBackgroundValue(back)
    edge_ob.SetForegroundValue(fore)
    edge_ob.SetMaximumError(maximErr)
    edge_ob.SetVariance(variance)

    print("Paramaters")
    print(edge_ob)

    print("Executing")
    out = edge_ob.Execute(img_f)

    print("saving")
    sitk.WriteImage(out, os.path.join(output, "zerocrossing.nrrd"))
    print("***Finished***")
    return (out, edge_ob)


def active_contour_level_set(seed, img, output, name="active_contour_ls.nrrd",
                             max_error=0.02, iterations=10000, curvature=1, propagation=1, advection=1):
    '''
    Harwell implementation of ITK GeodesicActiveContourLevelSetImageFilter
    @param seed = initial level set/seed (create using SignedMaurerDistanceMap)
    @param img = edge potential map
    @param output =  output folder
    @param GeodesicActiveContourLevelSetImageFilter arguments
    @return out = the segmented image
    @return ls_ob = the ls object
    '''
    print("***Harwell implementation of ITK GeodesicActiveContourLevelSetImageFilter***")
    print("setting up object")
    ls_ob = sitk.GeodesicActiveContourLevelSetImageFilter()
    ls_ob.SetMaximumRMSError(max_error)
    ls_ob.SetAdvectionScaling(advection)
    ls_ob.SetNumberOfIterations(iterations)
    ls_ob.SetCurvatureScaling(curvature)
    ls_ob.SetPropagationScaling(propagation)
    #lsFilter.ReverseExpansionDirectionOn()

    print("Paramaters:")
    print(ls_ob)

    print("Executing")
    out = ls_ob.Execute(seed, img)

    print("saving")
    sitk.WriteImage(out, os.path.join(output, name))
    print("***Finsihed***")
    return (out, ls_ob)


def threshold_level_set(seed, img, output, lower_threshold=200, upper_threshold=20000,
                        name="threshold_contour_ls.nrrd", max_error=0.02, iterations=10000,
                        curvature=1, propagation=1, advection=1):
    '''
    Harwell implementation of ITK ThresholdSegmentationLevelSetImageFilter
    @param seed = initial level set/seed (create using SignedMaurerDistanceMap)
    @param img = the image which is to be segmented
    @param output =  output folder
    @param ThresholdSegmentationLevelSetImageFilterarguments
    @return out = the segmented image
    @return ls_ob = the ls object
    '''

    print("***Harwell implementation of ITK ThresholdSegmentationLevelSetImageFilter***")
    print("setting up object")
    ls_ob = sitk.ThresholdSegmentationLevelSetImageFilter()
    ls_ob.SetLowerThreshold(lower_threshold)
    ls_ob.SetUpperThreshold(upper_threshold)
    ls_ob.SetMaximumRMSError(max_error)
    ls_ob.SetNumberOfIterations(iterations)
    ls_ob.SetCurvatureScaling(curvature)
    ls_ob.SetPropagationScaling(propagation)
    #lsFilter.ReverseExpansionDirectionOn()

    print("Paramaters:")
    print(ls_ob)

    print("Executing")
    out = ls_ob.Execute(seed, img)

    print("saving")
    sitk.WriteImage(out, os.path.join(output, name))
    print("***Finsihed***")
    return (out, ls_ob)


# Watershed seg
def ws_seg(img, output, sig=1.25, ws_level=0.9, force_background=False, watershed_line=True, fully_connected=False):
    # Get magnitude gradient
    print("Getting gaussian gradient magnitude")
    seg = sitk.GradientMagnitudeRecursiveGaussian(img, sigma=sig)
    sitk.WriteImage(seg, os.path.join(output, "gradient_4_ws.nrrd"))

    #sitk.WriteImage(seg,os.path.join(outputFolder,filter+"_blur.nrrd"))
    # Perform the watershed alg
    print("Doing watershed")
    seg = sitk.MorphologicalWatershed(seg, level=ws_level, markWatershedLine=watershed_line,
                                      fullyConnected=fully_connected)
    sitk.WriteImage(seg, os.path.join(output, "ws" + str(ws_level) + ".nrrd"))

    seg = sitk.RelabelComponent(seg)
    sitk.WriteImage(seg, os.path.join(output, "seg_componets_ws" + str(ws_level) + ".nrrd"))

    print("overlay")
    tif = os.path.join(output, "overlay_ws" + str(ws_level) + ".tif")
    if os.path.exists(tif):
        os.remove(tif)
    overlay = sitk.LabelOverlay(img, seg)
    sitk.WriteImage(overlay, os.path.join(output, "overlay_ws" + str(ws_level) + ".tif"))

    print("finished")
    return seg


# def watershed_mask_markers(img,output, sig = 1.25, ws_level = 0.9,force_background = False):
#     # Get magnitude gradient
#     print("Getting gaussian gradient magnitude"
#     seg = sitk.GradientMagnitudeRecursiveGaussian(img, sigma=sig)
#     sitk.WriteImage(seg,os.path.join(output,"gradient_4_ws.nrrd"))
#     size = img.GetSize()
#     print(size
#     pt = [size[0]/2, size[1]/2, size[2]/2]
#     print(pt
#
#     idx = img.TransformPhysicalPointToIndex(pt)
#
#     min_img = sitk.RegionalMinima(seg, backgroundValue=0, foregroundValue=1.0, fullyConnected=False, flatIsMinima=True)
#     sitk.WriteImage(min_img,os.path.join(output,"min_image.nrrd"))
#     marker_img = sitk.ConnectedComponent(min_img, fullyConnected=False)
#     sitk.WriteImage(seg,os.path.join(marker_img,"marker_image.nrrd"))
#
# #     marker_img[idx] = 2
#     #sitk.WriteImage(seg,os.path.join(outputFolder,filter+"_blur.nrrd"))
#     # Perform the watershed alg
# #     print("Doing watershed"
# #     ws = sitk.MorphologicalWatershedFromMarkers(seg, marker_img, markWatershedLine=True, fullyConnected=False)
# #     sitk.WriteImage(seg,os.path.join(output,"ws.nrrd"))
# #
# #     print("performing clean up"
# #     seg = clean_up(img,seg,output,"ws", force_background_clean = force_background)
# #
# #     print("finished"
# #     return seg


def edge_potential_map(img, output):
    # DOES NOT WORK!!!!!!!!!!!!!!!!!!!!!!!!
    print("***Get edge potential map DOES NOT WORK!***")
    print("blur")

    blur = sitk.GradientRecursiveGaussian(img, sigma=0.01)
    sitk.WriteImage(blur, os.path.join(output, "edge_potential_blur.nrrd"))

    edge = sitk.EdgePotential(blur)
    sitk.WriteImage(edge, os.path.join(output, "edge_potential_1.nrrd"))

    sig_ob = sitk.SigmoidImageFilter()

    #sig_ob.SetOutputMinimum( 0.0 );
    #sig_ob.SetOutputMaximum( 1.0 );
    #sig_ob.SetAlpha( -0.4 );
    #sig_ob.SetBeta( 2.5 );
    map = sig_ob.Execute(gradient)
    sitk.WriteImage(map, os.path.join(output, "edge_potential_map.nrrd"))

    return (map, sig_ob)
