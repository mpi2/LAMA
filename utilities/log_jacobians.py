#!/usr/bin/env python
import os
import SimpleITK as sitk
import numpy as np

def log_jacobians(indir, outdir):
    for path in os.listdir(indir):
        jac_path = os.path.join(indir, path)
        try:
            img = sitk.ReadImage(jac_path)
        except RuntimeError:
            continue
        arr = sitk.GetArrayFromImage(img)
        log_array = np.log(arr)
        out_img = sitk.GetImageFromArray(log_array)
        out_path = os.path.join(outdir, 'logjac_' + path)
        sitk.WriteImage(out_img, out_path)


if __name__ == '__main__':
    import sys
    indir = sys.argv[1]
    outdir = sys.argv[2]
    log_jacobians(indir, outdir)
