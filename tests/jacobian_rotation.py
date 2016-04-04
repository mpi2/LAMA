#!/usr/bin/env python

import math
import numpy as np
import SimpleITK as sitk
from scipy.linalg import sqrtm
import sys
sys.path.append('../utilities')
import transformations as trans


def run(invol):
    img = sitk.ReadImage(invol)
    arr = sitk.GetArrayFromImage(img)

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
                #scaling = s[0][:3]
                print r1
    # result_arr = np.array(result).reshape(arr.shape)



if __name__ == '__main__':
    import sys
    invol = sys.argv[1]
    run(invol)