"""
Test out whther we can do post processing on label segmetnations to improve lateral ventricle accuracy
"""

from skimage import measure
import SimpleITK as sitk
import numpy as np
from pathlib import Path
from typing import List


def write(obj, path):
    try:
        sitk.WriteImage(object, str(path))
    except NotImplementedError:
        sitk.WriteImage(sitk.GetImageFromArray(obj), str(path), True)


def overalap(a, a_label, b, b_label):
    a = np.copy(a)
    b = np.copy(b)
    a[a != a_label] = 0
    a[a == a_label] = 1
    b[b != b_label] = 0
    b[b == b_label] = 1
    return a[a == b].size


def run(image_to_segment: Path,
        initial_segmentation: Path):
    # Do right lateral ventricle first
    r_ventricle_label = 26
    r_surrounding_labels = [39, 29, 35]  # ventricular zone, thalamus, neopallial cortex and amygdala
    r_seg = segment_lateral_ventricles(image_to_segment, initial_segmentation, r_surrounding_labels, r_ventricle_label)

    # And the left
    l_ventricle_label = 25
    l_surrounding_labels = [40, 30, 36]
    l_seg = segment_lateral_ventricles(image_to_segment, initial_segmentation, l_surrounding_labels, l_ventricle_label)

    # merge the two segmentations
    seg = r_seg
    seg[l_seg == l_ventricle_label] = l_ventricle_label

    return seg


def segment_lateral_ventricles(image_to_segment_path: Path,
                               initial_segmentation_path: Path,
                               surrounding_labels: List[int],
                               target_label: int,
                               outpath=None) -> np.ndarray:
    """
    This is a lama segmetation add on function.
    It currently work well for the lateral ventricles in in E15.5 micro-CT images

    Parameters
    ----------
    image_to_segment_path:
        The image to be segmented
    initial_segmentation_path
        The inital LAMA segmentation
    surrounding_labels
        a list of label numbers that surround the target label
    target_label
        The label we are trying to improve

    Returns
    -------
    A new segmentation with the same dimensions as 'inital_segmentation' containing only the newly segmented
    'target_label'

    """

    image_to_segment = sitk.GetArrayFromImage(sitk.ReadImage(str(image_to_segment_path)))
    initial_segmentation = sitk.GetArrayFromImage(sitk.ReadImage(str(initial_segmentation_path)))

    # Remove all label that are not of interest here
    initial_segmentation[~np.isin(initial_segmentation, surrounding_labels + [target_label])] = 0

    # Set all the surrouding labels to 1
    initial_segmentation[np.isin(initial_segmentation, surrounding_labels)] = 1

    # write(lab, outdir / 'test_surroud.nrrd')

    props = measure.regionprops(initial_segmentation)

    for x in props:
        if x.label == 1:

            # Get the ROI defined by the surrounding labels
            b = x.bbox
            image_to_segment *= -1
            image_roi = image_to_segment[b[0]:b[3], b[1]: b[4], b[2]: b[5]]
            label_roi = initial_segmentation[b[0]:b[3], b[1]: b[4], b[2]: b[5]]

            # Do otsu thresholding on the roi
            ventr_area_itk = sitk.GetImageFromArray(image_roi)
            thresh = sitk.OtsuThreshold(ventr_area_itk)
            thresh_arr = sitk.GetArrayFromImage(thresh)
            thresh_arr = np.invert(thresh_arr)

            # Get connected component labels from the thresholding
            threshold_labels = measure.label(thresh_arr)
            new_segmentation = initial_segmentation.copy().astype(np.short)
            new_segmentation[:] = 0

            # Insert the thresgold label condidates back into the label map
            new_segmentation[b[0]:b[3], b[1]: b[4], b[2]: b[5]] = threshold_labels
            # write(new_segmentation, outdir / 'relabe_test.nrrd')

            # Find the threshold label with the largest overlap with the target label
            largest_overlap = 0
            largest_overlap_label = 0
            thresh_props = [x for x in measure.regionprops(threshold_labels) if x.area > 500]

            for pr in thresh_props:
                ol = overalap(label_roi, target_label, threshold_labels, pr.label)
                if ol > largest_overlap:
                    largest_overlap = ol
                    largest_overlap_label = pr.label

            new_segmentation[new_segmentation != largest_overlap_label] = 0
            new_segmentation[new_segmentation == largest_overlap_label] = target_label

            if outpath:
                write(new_segmentation.astype(np.uint8), outpath)

            return new_segmentation.astype(np.uint8)
