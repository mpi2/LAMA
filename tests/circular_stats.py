#!/usr/bin/env python


import numpy as np
import SimpleITK as sitk
import pycircstat as pc
import os
from os.path import join, basename, splitext


def run(indir, outdir):
    for name in os.listdir(indir):
        in_path = join(indir, name)
        img = sitk.ReadImage(in_path)
        arr = sitk.GetArrayFromImage(img)
        out_path = join(outdir, name)
        get_angles(arr, out_path)

def get_angles(arr, out_path):

    # Get the product of each jacobain matrix and its transform
    angles = []
    a = arr
    for b in a:
        for j in b:
           # calculate the angle between

            for vector in j:
                angle = np.rad2deg(np.arctan2(vector[0], vector[1]))
                angles.append(angle)

    result_arr = np.array(angles).reshape(arr.shape[0:3])
    result_img = sitk.GetImageFromArray(result_arr)
    sitk.WriteImage(result_img, out_path)



if __name__ == '__main__':
    import sys
    indir = sys.argv[1]#
    outdir = sys.argv[2]
    run(indir, outdir)



