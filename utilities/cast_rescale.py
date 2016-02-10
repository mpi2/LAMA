#!/usr/bin/env python

"""
Just for casting from 8 bit to 16bit with intesities scaled to max
"""

import SimpleITK as sitk
import numpy as np
from os.path import join


def run(img_path, out_path):
    img = sitk.ReadImage(img_path)
    cast = sitk.Cast(img, sitk.sitkUInt16)
    res = sitk.RescaleIntensity(cast, 0, 65536)
    sitk.WriteImage(res, out_path)


if __name__ == '__main__':
    import sys
    img_path = sys.argv[1]
    out_path = sys.argv[2]
    run(img_path, out_path)
