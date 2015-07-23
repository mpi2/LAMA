#!/usr/bin/python


"""
When finished, This module with take a labelmap and an 3D image. Then either dilate or shrink the region specified
 by the label.

 Only acepts one label at a time
"""
import SimpleITK as sitk
import numpy as np
import subprocess
from os.path import join
import os

ELASTIX_PATH = 'elxParamTemp.txt'

PARAM = """
(AutomaticScalesEstimation  "true")
(NumberOfHistogramBins  32)
(MovingImageDimension  3)
(DefaultPixelValue  0)
(UseDirectionCosines  "true")
(FixedImageDimension  3)
(NewSamplesEveryIteration  "true")
(FixedImagePyramid  "FixedSmoothingImagePyramid")
(FinalBSplineInterpolationOrder  3)
(Resampler  "DefaultResampler")
(WriteResultImage  "true")
(CompressResultImage  "true")
(ResultImageFormat  "nrrd")
(ImageSampler  "Random")
(WriteTransformParametersEachIteration  "false")
(FixedInternalImagePixelType  "short")
(MovingImagePyramid  "MovingSmoothingImagePyramid")
(Interpolator  "BSplineInterpolator")
(ResultImagePixelType  "float")
(AutomaticTransformInitialization  "true")
(MovingInternalImagePixelType  "short")
(HowToCombineTransforms  "Compose")
(ResampleInterpolator  "FinalBSplineInterpolator")
(NumberOfSpatialSamples  100000)
(Optimizer  "AdaptiveStochasticGradientDescent")
(FinalGridSpacingInVoxels  8)
(MaximumNumberOfIterations  100)
(ASGDParameterEstimationMethod  "DisplacementDistribution")
(Metric  "AdvancedNormalizedCorrelation")
(Registration  "MultiResolutionRegistration")
(Transform  "BSplineTransform")
(AutomaticParameterEstimation  "true")
(MaximumStepLength 1.0)
(NumberOfResolutions  1)
"""


RADIUS = 5
LABEL_VALUE = 200
FIXED = 'resized.nnrd'
MOVING = 'normal.nrrd'

def run(label, image):
    original, resized = resize_label(label, image)
    sitk.WriteImage(original, MOVING)
    sitk.WriteImage((resized, FIXED))

    reg_dir = 'out'
    os.mkdir(reg_dir)
    register_label_images(FIXED, MOVING, reg_dir)

    tform_param = join(reg_dir, 'TransformParameters.0.txt')


    warp_out = 'Warped'
    os.mkdir(warp_out)
    warp_organ(image, tform_param)



def resize_label(binary_label, image):

    label_img = sitk.ReadImage(binary_label)

    resized_label = sitk.BinaryDilate(label_img, RADIUS)

    img = sitk.ReadImage(image)
    orig_arr = sitk.GetArrayFromImage(img)
    resized_array = np.copy(orig_arr)
    orig_label = sitk.GetArrayFromImage(sitk.ReadImage(binary_label))
    resized_label = sitk.GetArrayFromImage(resized_label)

    orig_arr[orig_label == 1] = LABEL_VALUE
    resized_array[resized_label == 1] = LABEL_VALUE

    orig = sitk.GetImageFromArray(orig_arr)
    res = sitk.GetImageFromArray(resized_array)

    # Set the label in the iames to a number
    return orig, res

def register_label_images(original, enlarged, outdir):

        with open(ELASTIX_PATH, 'w') as fh:
            fh.write(PARAM)

        cmd = ['elastix',
               '-f', enlarged,
               '-m', original,
               '-out', outdir,
               '-p', ELASTIX_PATH
              ]

        subprocess.check_output(cmd)


def warp_organ(image, tform_file, out):

        cmd = ['transformix',
               '-in', image,
               '-tp', tform_file,
               '-out', out,
               ]

        subprocess.check_output(cmd)









if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Resize organs")
    parser.add_argument('-i', dest='image', required=True)
    parser.add_argument('-l', dest='label', required=True)
    args = parser.parse_args()

    run(args.image, args.label)
