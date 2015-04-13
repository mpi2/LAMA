#!/usr/bin/env python


import SimpleITK as sitk
import numpy as np


MIN = 0
MAX = 200


def run(volpath, out):
    vol = sitk.ReadImage(volpath)
    arr = sitk.GetArrayFromImage(vol)
    curr_min = arr.min()
    curr_max = arr.max()

    #Do a linear transform to move zero to 100
    im_s1 = sitk.RescaleIntensity(vol, curr_min + 100, curr_max + 100)

    # Recale so range goes from 0-200
    im_s2 = sitk.RescaleIntensity(im_s1, MIN, MAX)

    sitk.WriteImage(im_s2, out)


if __name__ == '__main__':
    import sys
    vol = sys.argv[1]
    out = sys.argv[2]

    run(vol, out)
