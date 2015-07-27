#!/usr/bin/env python

"""
This is just a temp file to test out the GLCM code on some mutants.

If it works, incorpoate  into reg_stats.py. But make reg_stats a class and subclass for each stats method as it's
getting a bit out of hand
"""
from scipy import stats

import harwellimglib as hil
from utilities.glcm3d import Glcm
import numpy as np
import SimpleITK as sitk


def run(wt_dir, mut_dir, output_img, chunksize):
    wts = hil.GetFilePaths(wt_dir)
    muts = hil.GetFilePaths(mut_dir)

    shape = sitk.GetArrayFromImage(sitk.ReadImage(wts[0])).shape

    print 'getting wt glcms'
    wt_contrasts = []
    for wt in wts:
        glcm_maker = Glcm(wt, chunksize)
        wt_contrasts.append(glcm_maker.get_contrasts())

    print "getting mut glcms"
    mut_contrasts = []
    for mut in muts:
        glcm_maker = Glcm(mut, chunksize)
        mut_contrasts.append(glcm_maker.get_contrasts())

    print 'doing stats'
    print 'doing stats'
    wt_stacked = np.vstack(wt_contrasts)
    mut_stacked = np.vstack(mut_contrasts)

    raw_stats = stats.ttest_ind(wt_stacked, mut_stacked)

    # reform a 3D array from the stas and write the image
    out_array = np.zeros(shape)

    i = 0

    for z in range(0, shape[0] - chunksize, chunksize):
        print 'w', z
        for y in range(0, shape[1] - chunksize, chunksize):
            for x in range(0, shape[2] - chunksize, chunksize):
                score = raw_stats[0][i]
                prob = raw_stats[1][i]
                if prob < 0.05:
                    output_value = score
                else:
                    output_value = 0
                out_array[z: z + chunksize, y: y + chunksize, x: x + chunksize] = output_value
                i += 1

    out = sitk.GetImageFromArray(out_array)
    sitk.WriteImage(out, output_img)





if __name__ == '__main__':
    import sys
    wt_dir = sys.argv[1]
    mut_dir = sys.argv[2]
    output_img = sys.argv[3]
    chunksize = sys.argv[4]
    run(wt_dir, mut_dir, output_img, chunksize)
