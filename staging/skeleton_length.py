#!/usr/bin/env python


from scipy.ndimage.measurements import center_of_mass
from os.path import join
import os
import SimpleITK as sitk
import numpy as np
from scipy.spatial.distance import cdist


def run(in_dir, verbose=False):
    lengths = {}
    for path, subdirs, files in os.walk(in_dir):
        for name in files:
            if not name.endswith('nrrd'):
                continue
            im_path = join(path, name)
            img = sitk.ReadImage(im_path)
            arr = sitk.GetArrayFromImage(img)
            dist = skeletonize(arr)
            # Bodge: Remove any se_ prefixes from the inverted segmentation
            name = name.strip('seg_')
            if verbose:
                print(("{},{}".format(name, dist)))
            lengths[name] = dist
    return lengths


def skeletonize(arr):
    out_arr = np.zeros_like(arr)
    points = []
    for z, slice_ in enumerate(arr):
        if np.any(slice_):
            center = center_of_mass(slice_)
            try:
                out_arr[z, int(center[0]), int(center[1])] = 1
            except IndexError:
                print('ogh no')
            points.append((z, center[0], center[1]))
    # created two padded arrays so they are shifted relative to each other
    points_1 = points[:]
    points_1.insert(0, points[0])
    points_2 = points[:]
    points_2.append(points[-1])
    # Pairs that we want are on the diagonal
    distances = cdist(points_1, points_2)
    adjacent_distances = np.diag(distances)
    total_distance = np.sum(adjacent_distances)
    return total_distance


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("?")
    parser.add_argument('-i', '--indir', dest='indir', help='Input path.', required=True)
    args = parser.parse_args()

    run(args.indir, True)
