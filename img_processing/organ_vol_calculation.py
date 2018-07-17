#! /usr/bin/env python3

"""
See lab book dated 5th/6th June 18

Get organ volumes from a bunch or inverted label maps
"""

import os
from os.path import join, split
import addict
import SimpleITK as sitk
import pandas as pd
from common import get_file_paths

def normalised_label_sizes(label_dir, mask_dir, outpath):
    """
    Given a directory of labelmaps and whole embryo masks, generate a csv file containing organ volumes normalised to
    mask size

    Parameters
    ----------
    label_dir: str
        directory containing (inverted) label maps - can be in subdirectories
    mask_dir: str
        directory containing (inverted) masks - can be in subdirectories
    outpath: str
        path to save generated csv

    """

    label_df = _get_label_sizes(get_file_paths(label_dir))

    mask_df = _get_label_sizes(get_file_paths(mask_dir))

    normed_label_volumes =  label_df.divide(mask_df[1], axis=0)

    normed_label_volumes.to_csv(outpath)


def _get_label_sizes(paths):
    """
    Get the organ volumes for a bunch of of speciemns
    Parameters
    ----------
    label_paths: list
        paths to labelmap volumes

    Returns
    -------
    pandas dataframe:
        columns: label (organ)
        rows: specimen ids
    """

    label_volumes = addict.Dict()
    to_do = len(paths)
    n = 1

    for label_path in paths:

        print("{} of {}".format(n, to_do))
        n += 1
        # Get the name of the volume
        volname = os.path.split(split(label_path)[0])[1]
        labelmap = sitk.ReadImage(label_path)
        arr = sitk.GetArrayFromImage(labelmap)
        max_label = int(arr.max())


        lsf = sitk.LabelStatisticsImageFilter()
        labelmap = sitk.Cast(labelmap, sitk.sitkUInt16)
        lsf.Execute(labelmap, labelmap)
        # max_label = lsf.GetMaximum()
        for i in range(1, max_label + 1):
            voxel_count = lsf.GetCount(i)
            label_volumes[volname][i] = voxel_count
    return pd.DataFrame(label_volumes.to_dict()).T


if __name__ == '__main__':
    pass

    # home = os.path.join(os.path.expanduser('~'), 'bit')
    #
    # invert_labels_dir = '/home/neil/bit/LAMA_results/E14.5/paper_runs/mutant_runs/280618_analysed_lines/ATG13/output/inverted_labels/similarity'
    # inverted_mask_dir = '/home/neil/bit/LAMA_results/E14.5/paper_runs/mutant_runs/280618_analysed_lines/ATG13/output/inverted_stats_masks/similarity'
    # outdir = '/home/neil/bit/LAMA_results/E14.5/paper_runs/mutant_runs/280618_analysed_lines/ATG13/output/inverted_labels'
    #
    #
    # dirs = os.listdir(invert_labels_dir)
    #
    # label_df = get_label_vols(get_file_paths(invert_labels_dir)).T
    # label_df.to_csv(join(outdir, '050718_raw_label_vols.csv'))
    #
    # mask_df = get_label_vols(get_file_paths(inverted_mask_dir)).T
    # mask_df.to_csv(join(outdir, '050718_mask_vols.csv'))
    #
    # normed_label_volumes =  label_df.divide(mask_df[1], axis=0)
    # normed_outpath = join(outdir, 'organ_volumes_normed_to_mask.csv')
    # normed_label_volumes.to_csv(normed_outpath)
