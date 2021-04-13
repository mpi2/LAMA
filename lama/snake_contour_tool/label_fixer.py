"""
fixes poor labels specified within a config file (e.g. fixer_config.toml):

config file should contain the rigidly alligned scan and the inverted labels in separate directories
specified by rigid_img_path and inverted_labels_path

parameters for config file (each column is a separate run):
lab_of_int = label to be corrected

kmean = number of clusters generated to seperate organs of interest

num_int = number of interations to grow/shrink label to the cluster

prop_scaling = scaling factor (negative to remove overlabelling, positive for labelling)

example run:

python3 label_fixer.py -c fixer_config.toml
"""
from pathlib import Path
import nrrd
import SimpleITK as sitk
import toml

from lama import common
from lama.snake_contour_tool import geodesic_active_segmentation
import numpy as np
import argparse


def get_config(config_file):
    """Gets a csv file containing the order of labels to fisx and their numb_int and kmean_vals """
    config = toml.load(config_file)
    return config


def cluster(rigid_img, i):
    """Do the kmeans clustering, i is the only parameter for this kmeans function"""
    k = sitk.ScalarImageKmeans(rigid_img, [i] * i)

    out_arr = sitk.GetArrayFromImage(k)
    return out_arr


def fix_label(rigid_img_path, inverted_labels_path, label_of_interest, kmean, num_iterations, prop_scaling):
    """Uses the active contour tool to automatically segment a poorly inverted label to a cluster image"""

    rigid, r_head = nrrd.read(rigid_img_path)
    inv_labels, inv_head = nrrd.read(inverted_labels_path)

    # # Get the bounding box of the inverted label - # Kyle - you don't need an ROI to fix the segmentation.
    #                                               # just get the label of interest
    single_label, _ = nrrd.read(inverted_labels_path)
    single_label[single_label != label_of_interest] = 0

    # Perform the snake contouring to get a corrected label
    rigid_img = sitk.GetImageFromArray(rigid)
    out_arr = cluster(rigid_img, kmean)

    corrected_label = geodesic_active_segmentation.snake(single_label, out_arr, label_of_interest, num_iterations,
                                                         prop_scaling)

    # Replace the old label with the corrected label
    corrected_array = sitk.GetArrayFromImage(corrected_label)

    # set over-labelled vals in original label to 0
    overlabelled_vals = np.logical_and(inv_labels == label_of_interest, corrected_array != 1)
    inv_labels[overlabelled_vals == True] = 0

    # set non-labelled values to the label of interest and write the file
    inv_labels[corrected_array == 1] = label_of_interest
    nrrd.write(str(inverted_labels_path), inv_labels, header=inv_head)


def main():
    parser = argparse.ArgumentParser("fix specified labels using the active contour tool")
    parser.add_argument('-c', '--config', dest='config_file', help='config file for label_fixer')

    args = parser.parse_args()

    # get configuration by reading the toml file
    config = get_config(args.config_file)

    # Set the rigid img and inverted label paths from the config
    rigid_img_path = config["rigid_img_path"]
    inverted_labels_path = config["inverted_labels_path"]

    # there should be multiple rigid images and labels, therefore get all filepaths
    rigid_paths = common.get_file_paths(Path(rigid_img_path))
    label_paths = common.get_file_paths(Path(inverted_labels_path))

    # for each volume, perform multiple active contours in the order and with the parameters specified within the config
    for i in range(len(rigid_paths)):
        parameters = config['active_contour_params']
        for (lab_of_int, kmean, num_int, prop_scaling) in zip(parameters['lab_of_int'], parameters['kmean'],
                                                              parameters['num_int'],
                                                              parameters['prop_scaling']):
            fix_label(Path(rigid_paths[i]), Path(label_paths[i]), lab_of_int, kmean, num_int, prop_scaling)


if __name__ == '__main__':
    main()
