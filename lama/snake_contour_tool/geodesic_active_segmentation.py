#Kyle Drover - Geodesic Active Contour Segmentation based on existing label

from pathlib import Path
import SimpleITK as sitk
import nrrd
from scipy import ndimage
import numpy as np

def snake(label, cluster, label_of_interest, N, propScaling):

    # Set the distance
    distance = sitk.SignedMaurerDistanceMapImageFilter()
    distance.InsideIsPositiveOff()
    distance.UseImageSpacingOn()

    # set the seed from the label and make the initial image from the seed
    seedImage = sitk.GetImageFromArray(label.astype(np.int16), isVector=False)

    initialImage = sitk.BinaryThreshold(distance.Execute(seedImage), -1000, 10)
    initialImage = sitk.Cast(initialImage, sitk.sitkFloat32) * -1 + 0.5

    # Setting up the feature image

    # determine which cluster has the maximum overlay with the label.
    c_scores = []
    for k in range(np.amax(cluster)):
        overlay = np.logical_and(label == label_of_interest, cluster == k)
        c_scores.append(np.sum(label[overlay == True]))


    c_val = np.argmax(c_scores)

    # convert desired cluster into a binary image
    cluster[cluster == c_val] = 8


    #cluster[cluster != label_of_interest] = 0
    #cluster[cluster == label_of_interest] = 10

    # turn cluster into sitk image and blur it
    featureImage = sitk.GetImageFromArray(cluster.astype(np.int16), isVector=False)

    gradientMagnitude = sitk.GradientMagnitudeRecursiveGaussianImageFilter()# feature image needs blurring for some reason
    gradientMagnitude.SetSigma(0.05)# sigma is the blurring factor

    featureImage = sitk.BoundedReciprocal(gradientMagnitude.Execute(featureImage))
    # featureImage = sitk.InvertIntensity(featureImage, 1)
    featureImage = sitk.Cast(featureImage, sitk.sitkFloat32)

    test_feat_arr= sitk.GetArrayFromImage(featureImage)
    nrrd.write(
        'C:/Users/u5823099/Anaconda3/Lib/site-packages/lama/tests/test_data/reg_fixer_data/inverted_masks/2086980_feat_arr.nrrd',
        test_feat_arr)

    # sitk.WriteImage(featureImage, 'C:/Users/u5823099/Anaconda3/Lib/site-packages/lama/tests/test_data/reg_fixer_data/inverted_masks/2086980_feature.nrrd')
    # perform the active contour
    geodesicActiveContour = sitk.GeodesicActiveContourLevelSetImageFilter()
    geodesicActiveContour.SetPropagationScaling(propScaling)
    geodesicActiveContour.SetCurvatureScaling(0.15)
    geodesicActiveContour.SetAdvectionScaling(0.15)
    geodesicActiveContour.SetMaximumRMSError(0.01)
    geodesicActiveContour.SetNumberOfIterations(N)

    levelset = geodesicActiveContour.Execute(initialImage, featureImage)

    #convert the output of the active contour (i.e. the levelset) to a binary image
    bi_image = sitk.BinaryThreshold(levelset, -1000, 0)

    return bi_image