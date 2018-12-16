#!/usr/bin/env python

import SimpleITK as sitk
import os
from lama import common


def batch_hm(dir_, target_path, out_dir, rescale=False):
    for in_path in common.get_file_paths(dir_):

        matched = hm(in_path, target_path)
        base = os.path.basename(in_path)
        out_path = os.path.join(out_dir, base)
        sitk.WriteImage(matched, out_path)


def hm(input_path, target_path):

    input_img = sitk.ReadImage(input_path)
    target_img = sitk.ReadImage(target_path)
    matched = sitk.HistogramMatching(input_img, target_img)
    return matched


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser("Histogram match two images")
    parser.add_argument('-i', '--input', dest='input_dir', help='Image dir which will have intensity matched', required=True)
    parser.add_argument('-t', '--target', dest='target', help='Image to get target histogram from', required=True)
    parser.add_argument('-o', '--out', dest='out_dir', help='Image to get target histogram from', required=True)
    parser.add_argument('-r', '--rescale_first', dest='rescale', help='Rescale to 0-255 before matching', default=False,
                        action='store_true')

    args = parser.parse_args()

    batch_hm(args.input_dir, args.target, args.out_dir, args.rescale)


