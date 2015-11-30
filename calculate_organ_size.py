#!/usr/bin/env python

import csv
import common
from os.path import join, split, basename
import SimpleITK as sitk
import numpy as np


def calculate_volumes(labels_dir, label_names, outfile, voxel_size):
    """
    Calculate organ volumes and save results in a csv file

    Parameters
    ----------
    labels_dir: str
        path to dir containing label maps. Volumes can exist in subdirectories.
    label_names: str
        path to csv file containing label names:
        eg:
            liver
            thyroid
            left kidney
    outfile: str
        path to save calculated volumes as csv file
    voxel_size: float
        size of voxels
    """
    voxel_size = float(voxel_size)

    header = [' ']
    if label_names:
        with open(label_names, 'r') as headerfile:
            reader = csv.reader(headerfile)
            header.extend(reader.next())

    with open(outfile, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(header)
        labelmaps = common.GetFilePaths(labels_dir)
        for label_path in labelmaps:
            labelmap = sitk.ReadImage(label_path)
            label_array = sitk.GetArrayFromImage(labelmap)
            hist = np.histogram(label_array, bins=label_array.max() + 1)[0] # How works this?
            num_voxels = list(hist)
            vols = [n * voxel_size for n in num_voxels]  # to get um^3 ?

            # Get the name of the volume
            volname = split(split(label_path)[0])[1]
            # Remove first value. as this is the background and replace with vol name
            vols[0] = basename(volname)

            writer.writerow(vols)

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser("Calculate organ volumes")
    parser.add_argument('-l', '--labels_dir', dest='labels_dir', help='Folder containing inverted labels (can be in sub dirs)', required=True)
    parser.add_argument('-n', '--names', dest='label_names', help='csv file with label names, one per row', required=True)
    parser.add_argument('-o', '--out', dest='out_file', help='Outfile csv', required=True)
    parser.add_argument('-v', '--vsize', dest='voxel_size', help='Voxel size of specimens', default=False)
    args, _ = parser.parse_known_args()
    calculate_volumes(args.labels_dir, args.label_names, args.out_file, args.voxel_size)
