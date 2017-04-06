#!/usr/bin/env python

import SimpleITK as sitk
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
import common


def convert_16_bit_to_8bit(indir, outdir):

    paths = common.GetFilePaths(indir)

    for path in paths:
        img = sitk.ReadImage(path)
        arr = sitk.GetArrayFromImage(img)
        # omin = arr.min()
        # omax = arr.max()
        # base = os.path.basename(path)
        if arr.dtype in (np.uint8, np.int8):
            # print("skipping {}. Already 8bit".format(base))
            continue

        arr2 = arr/256
        print arr2.min(), arr2.max()
        # cmin = arr.min()
        # cmax = arr.max()
        # print "{}. original min-max {} {}. Converted min-max: {} {}".format(base, omin, omax, cmin, cmax)
        #arr = arr.astype(np.uint8)
        out_img = sitk.GetImageFromArray(arr2)
        basename = os.path.basename(path)
        outpath = os.path.join(outdir, basename)
        sitk.WriteImage(out_img, outpath, True)


def cast_and_rescale_to_8bit(indir, outdir):

    paths = common.GetFilePaths(indir)

    for path in paths:
        img = sitk.ReadImage(path)
        cast = sitk.Cast(sitk.RescaleIntensity(img), sitk.sitkUInt8)
        basename = os.path.basename(path)
        outpath = os.path.join(outdir, basename)
        sitk.WriteImage(cast, outpath, True)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser("Rescale 16 bit images to 8bit")
    parser.add_argument('-i', dest='indir', help='dir with vols to convert', required=True)
    parser.add_argument('-o', dest='outdir', help='dir to put vols in', required=True)
    args = parser.parse_args()
    convert_16_bit_to_8bit(args.indir, args.outdir)