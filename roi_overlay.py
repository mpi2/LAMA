#!/usr/bin/env python


from os.path import splitext, basename, join
import SimpleITK as sitk
import numpy as np
from scipy.misc import imsave
import common
import logging

def make_normalization_roi_qc_images(img_dir, roi, out_dir):
    """
    Make some QC images showing the roi used for normalization overlaid on the registered, normalised images
    Parameters
    ----------
    norm_image_folder: str
        path to normalised images
    roi: list
        [roi starts, roi_ends] z,yx
    """

    file_paths = common.GetFilePaths(img_dir)
    if not file_paths or len(file_paths) < 1:
        return

    roi_starts, roi_ends = roi
    for img_path in file_paths:

        img = sitk.ReadImage(img_path)
        cast_img = sitk.Cast(sitk.RescaleIntensity(img), sitk.sitkUInt8)
        arr = sitk.GetArrayFromImage(cast_img)

        try:
            sag_slice_index = roi_starts[2] + ((roi_ends[2] - (roi_starts[2])) /2)
            cor_slice_index = roi_starts[1] + ((roi_ends[1] - (roi_starts[1])) / 2)
            ax_slice_index = roi_starts[0] + ((roi_ends[0] - (roi_starts[0])) / 2)

            sag_slice = arr[:, :, sag_slice_index]
            cor_slice = arr[:, cor_slice_index, :]
            ax_slice = arr[ax_slice_index, :, :]
        except IndexError:
            logging.warn("Cannot generate roi QC overlays. ROi is out of bounds")
            return
        roi_props = []
        roi_props.append([sag_slice, roi_starts[0:2], roi_ends[0:2], True])
        roi_props.append([cor_slice, [roi_starts[0], roi_starts[2]], [roi_ends[0], roi_ends[2]], True])
        roi_props.append([ax_slice,  [roi_starts[1], roi_starts[2]], [roi_ends[1], roi_ends[2]], False])
        # Draw roi on the slice
        widths = []
        heights = []
        images = []
        for slice_, roi_starts, roi_ends, do_flip in roi_props:
            widths.append(slice_.shape[1])
            heights.append(slice_.shape[0])
            yellow_indices = bounding_box_indices(roi_starts, roi_ends)
            rgb_arr = grey_to_rgb(slice_)
            for index in yellow_indices:
                rgb_arr[index[0], index[1]] = [255, 255, 0]
            if do_flip:
                images.append(np.flipud(rgb_arr))
            else:
                images.append(rgb_arr)

        # Vreate an image array the combined width of the threee images
        max_height = max(heights)
        total_width = sum(widths)
        out_img_arr = np.zeros(max_height * total_width * 3).reshape((max_height, total_width, 3))

        accumulated_width = 0

        for single_rgb_img in images:
            width = single_rgb_img.shape[1]
            height = single_rgb_img.shape[0]
            print out_img_arr[0: height, accumulated_width:width].shape
            print single_rgb_img.shape
            out_img_arr[0: height, accumulated_width: accumulated_width + width] = single_rgb_img
            accumulated_width += width
        base = splitext(basename(img_path))[0]
        out_path = join(out_dir, base + '.png')
        imsave(out_path, out_img_arr)


def bounding_box_indices(roi_starts, roi_ends):
    indices = []
    y1 = roi_starts[0]
    y2 = roi_ends[0]
    x1 = roi_starts[1]
    x2 = roi_ends[1]

    #left vertical
    for i in range(y1, y2):
        indices.append([i, x1])
    # right vertical
    for i in range(y1, y2):
        indices.append([i, x2])
    # top row
    for i in range(x1, x2):
        indices.append([y1, i])
    # bottom row
    for i in range(x1, x2):
        indices.append([y2, i])
    return indices



def grey_to_rgb(im):
    return np.asarray(np.dstack((im, im, im)), dtype=np.uint8)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("The MRC Harwell image registration pipeline")
    parser.add_argument('-i', dest='in_img', help='Config file (YAML format)', required=True)
    parser.add_argument('-o', dest='out_dir', help='Config file (YAML format)', required=True)
    args = parser.parse_args()

    roi = [[146, 109, 67], [156, 119, 77]]

    make_normalization_roi_qc_images(args.in_img, roi, args.out_dir)