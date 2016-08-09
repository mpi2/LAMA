#!/usr/bin/env python

"""
This module was created in order to extract 3D ROIs from full resolution CT/HREM images.
The input is a folder of z slices and an roi of interest.
It outputs a 3D nrrd image
"""


import common
import numpy as np
import SimpleITK as sitk
from PIL import Image


def write_roi_from_image_stack(indir, outpath, starts, ends, pattern=None):
    """

    Parameters
    ----------
    indir
    outdir
    starts: list
        the roi starting coords (xyz)
    ends: list
        the roi end coords (xyz)
    pattern

    Returns
    -------

    """
    x1, y1, z1 = [int(x) for x in starts]
    x2, y2, z2 = [int(x) for x in ends]
    depth = z2 - z1
    width = x2 - x1
    height = y2 - y1

    file_list = list(sorted(common.GetFilePaths(indir)))
    if pattern:
        file_list = [x for x in file_list if pattern in x]

    # SimpleITK reads in 2D bmps as 3D. So use PIL instead
    if file_list[0].lower().endswith('.bmp'):
        reader = pil_load
    else:
        reader = sitk_load

    out_array = np.zeros(depth * width * height).reshape((depth, height, width))

    z_slices = file_list[z1 + 1:z2]
    for i, path in enumerate(z_slices):
        slice_ = reader(path)
        roi_slice = slice_[y1:y2, x1:x2]
        out_array[i, :, :] = roi_slice

    roi_img = sitk.GetImageFromArray(out_array)
    sitk.WriteImage(roi_img, outpath)


def sitk_load(p):
    img = sitk.ReadImage(p)
    return sitk.GetArrayFromImage(img)


def pil_load(p):
    im = Image.open(p)
    return np.array(im)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser("The MRC Harwell phenotype detection tool")
    parser.add_argument('-in', dest='indir', help='dir with vols to convert', required=True)
    parser.add_argument('-out', dest='outdir', help='dir to put vols in', required=True)
    parser.add_argument('-s', dest='starts', help='start indices (x, y, z)', required=True, nargs=3, type=int)
    parser.add_argument('-e', dest='ends', help='end indices (x, y, z)', required=True, nargs=3, type=int)
    parser.add_argument('-p', dest='pattern', help="file name pattern eg '_rec'", required=False, type=str, default=None)
    args = parser.parse_args()
    write_roi_from_image_stack(args.indir, args.outdir, args.starts, args.ends, args.pattern)
