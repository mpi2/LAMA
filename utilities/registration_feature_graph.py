#! /usr/bin/env python

"""

"""

import SimpleITK as sitk
from SimpleITK import ReadImage, WriteImage, GetArrayFromImage, GetImageFromArray
import numpy as np
from os.path import join
import os
import sys
import matplotlib.pyplot as plt
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common
from scipy.spatial.distance import dice, jaccard



def run(small_path, large_path):
    # im = ReadImage(lama_data_path)
    # arr = GetArrayFromImage(im)

    # Try out gettiong features along z axis
    small_im = sitk.ReadImage(small_path)
    large_im = sitk.ReadImage(large_path)

    small = sitk.GetArrayFromImage(small_im)
    large = sitk.GetArrayFromImage(large_im)

    # small[small > 0] = 0
    # large[large > 0] = 0
    #
    # small[small < 0] = 1
    # large[large < 0] = 1

    small[small < 0] = 0
    large[large < 0] = 0

    small[small > 0] = 1
    large[large > 0] = 1

    shape = small.shape
    chunk_size = 8

    small_out = []
    large_out = []
    dice_result = 1 - dice(small.flatten().astype(np.uint8), large.flatten().astype(np.uint8))
    print dice_result

    # for x in range(0, shape[2] - chunk_size, chunk_size):
    #     for y in range(0, shape[1] - chunk_size, chunk_size):
    #         for z in range(0, shape[0] - chunk_size, chunk_size):
    #
    #             small_chunk = small[z: z + chunk_size, y: y + chunk_size, x: x + chunk_size]
    #             large_chunk = large[z: z + chunk_size, y: y + chunk_size, x: x + chunk_size]
    #
    #             if np.any(small_chunk):
    #                 small_out.append(1)
    #             else:
    #                 small_out.append(0)
    #             if np.any(large_chunk):
    #                 large_out.append(1)
    #             else:
    #                 large_out.append(0)
    #
    # dice_result = 1- dice(np.array(small_out), np.array(large_out))
    #print dice_result





if __name__ == '__main__':


    small = '/home/neil/sig/LAMA_results/E14.5/target_testing/090516_target_testing/comparison/RBMS1/small_rigid_flipped_rbms1.nrrd'
    large = '/home/neil/sig/LAMA_results/E14.5/target_testing/090516_target_testing/comparison/RBMS1/large_rigid_flipped__rbms1.nrrd'

    run(small, large)
