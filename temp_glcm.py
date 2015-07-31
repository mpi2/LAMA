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


def run(wt_dir, mut_dir, chunksize, mask):
    wts = hil.GetFilePaths(wt_dir)
    muts = hil.GetFilePaths(mut_dir)

    shape = sitk.GetArrayFromImage(sitk.ReadImage(wts[0])).shape

    print "getting mut glcms"
    mut_glcms = []
    for mut in muts:
        glcm_maker = Glcm(mut, chunksize, mask)
        mut_glcms.append(glcm_maker.get_glcms())
    np.save('mut_glcms_5px.npy', mut_glcms)

    print 'getting wt glcms'
    wt_glcms = []
    for wt in wts:
        glcm_maker = Glcm(wt, chunksize, mask)
        wt_glcms.append(glcm_maker.get_glcms())
    np.save('wt_glcms_5px.npy', wt_glcms)



    return
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
    chunksize = int(sys.argv[3])
    mask = sys.argv[4]
    run(wt_dir, mut_dir, chunksize, mask)
