#!/usr/bin/env python

import math
import numpy as np
import SimpleITK as sitk
from scipy.linalg import sqrtm
import sys
sys.path.append('..')
import transformations as trans
from os.path import join
import os



def run(indir, outdir):
    for name in os.listdir(indir):
        in_path = join(indir, name)
        img = sitk.ReadImage(in_path)
        arr = sitk.GetArrayFromImage(img)
        out_path = join(outdir, name)
        get_angles(arr, out_path)

def get_angles(arr, out_path):

    # Get the product of each jacobain matrix and its transform
    result = []
    for i, a in enumerate(arr):
        if i % 1000 != 0:
            continue
        for b in a:
            for j in b:
                j = j.reshape((3, 3))
                jT = j.T
                vP = np.dot(jT, j)
                s = sqrtm(vP)
                r = np.dot(np.linalg.inv(s), j)
                result.append(r)

                # Try to see if R^-1 == R.T
                #rT = r.T
                r1 = np.linalg.inv(r)
                #rD = np.linalg.det(r)
                sE = np.linalg.eig(s)[1][0]
                ddv =  np.dot(r, sE)
                row = [0,0,0]
                r = np.vstack((r, row))
                col = np.array([0, 0, 0, 1])
                r = np.vstack((r.T, col)).T
                scale = trans.decompose_matrix(r)
                angles = np.rad2deg(scale[2])
                angle = angles[0]
                result.append(angle)

    result_arr = np.array(result).reshape(arr.shape)
    sitk.WriteImage(result_arr, out_path)



if __name__ == '__main__':
    import sys
    indir = sys.argv[1]#
    outdir = sys.argv[2]
    run(indir, outdir)