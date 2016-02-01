#!/usr/bin/env python

import csv
import common
from os.path import join, split, basename
import SimpleITK as sitk
import logging


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

    # get the label names
    header = [' ']
    if label_names:
        with open(label_names, 'r') as name_reader:
            for line in name_reader:
                label_name = line.strip()
                header.append(label_name)

    with open(outfile, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(header)
        labelmaps = common.GetFilePaths(labels_dir)
        for label_path in labelmaps:
            # Get the name of the volume
            volname = split(split(label_path)[0])[1]
            speciemn_row = [volname]
            labelmap = sitk.ReadImage(label_path)
            lsf = sitk.LabelStatisticsImageFilter()
            lsf.Execute(labelmap, labelmap)
            num_labels = lsf.GetNumberOfLabels()
            if num_labels != len(header):
                logging.warn('The number of labels in inverted labelmap({}) does not match number of label names in {}'
                             .format(label_path, label_names))
            for i in range(1, num_labels):  # skip 0: unlabelled regions
                voxels = lsf.GetCount(i)
                speciemn_row.append(str(voxels))

            writer.writerow(speciemn_row)

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser("Calculate organ volumes")
    parser.add_argument('-l', '--labels_dir', dest='labels_dir', help='Folder containing inverted labels (can be in sub dirs)', required=True)
    parser.add_argument('-n', '--names', dest='label_names', help='csv file with label names, one per row', required=True)
    parser.add_argument('-o', '--out', dest='out_file', help='Outfile csv', required=True)
    parser.add_argument('-v', '--vsize', dest='voxel_size', help='Voxel size of specimens', default=False)
    args, _ = parser.parse_known_args()
    calculate_volumes(args.labels_dir, args.label_names, args.out_file, args.voxel_size)
