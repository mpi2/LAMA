"""Normalises the radiomics scans by the average intensity of a mask"""
from typing import Union

from pandas import Series, DataFrame

from lama.img_processing import normalise
from logzero import logger as logging
from lama import common
import os
import nrrd
from pathlib import Path

import numpy as np
import SimpleITK as sitk

from radiomics import featureextractor

import pandas as pd


# each scan in Ben's dataset will need its own mask
def get_images_from_masks(dir):
    img_list = []
    spec_name_list = []
    mask_list = []
    scan_paths = [spec_path for spec_path in common.get_file_paths(dir) if ('imgs' in str(spec_path))]
    mask_paths = [mask_path for mask_path in common.get_file_paths(dir) if ('labels' in str(mask_path))]

    # enumerate for indexing masks
    for i, img_path in enumerate(scan_paths):
        mask, m_h = nrrd.read(mask_paths[i])

        img, img_h = nrrd.read(img_path)
        # Only get values inside of the mask
        logging.info(f"Obtaining values from {img_path}")

        # s = ndimage.find_objects(mask)[0]

        # mask = mask[s[0].start:s[0].stop,
        #       s[1].start:s[1].stop,
        #       s[2].start:s[2].stop]
        # img = img[s[0].start:s[0].stop,
        #      s[1].start:s[1].stop,
        #      s[2].start:s[2].stop]

        img[(mask != 1) | (img < 0)] = 0

        spec_name_list.append(os.path.splitext(img_path.name)[0])
        # print(spec_name_list)
        img_list.append(img)
        mask_list.append(mask)
    return img_list, spec_name_list, mask_list


def pyr_calc_all_features(dir, normed: bool = False):
    # get either the normalised or original images
    scan_paths = [spec_path for spec_path in common.get_file_paths(dir) if ('normed' in str(spec_path))] if normed \
        else [spec_path for spec_path in common.get_file_paths(dir) if ('imgs' in str(spec_path))]

    tumour_paths = [spec_path for spec_path in common.get_file_paths(dir) if ('tumour_respaced' in str(spec_path))]

    # debugging - Thanks Neil
    scan_paths.sort()
    tumour_paths.sort()

    # Get the first order measurements
    full_orders = []

    for i, img_path in enumerate(scan_paths):
        logging.info(f"Calculating for {os.path.splitext(os.path.basename(img_path))[0]}")
        # cant use sitk.readimage due to header?
        img, img_head = nrrd.read(img_path)
        img = sitk.GetImageFromArray(img)
        mask, m_head = nrrd.read(tumour_paths[i])
        mask = sitk.GetImageFromArray(mask)

        # get all features and append to list
        extractor = featureextractor.RadiomicsFeatureExtractor()

        result = extractor.execute(img, mask)

        first_orders = pd.DataFrame.from_dict(result, orient='index',
                                              columns=[os.path.splitext(os.path.basename(img_path))[0]])
        full_orders.append(first_orders)

    return pd.concat(full_orders, axis=1)


def main():
    logging.info("Calculating Original First Order Features")
    _dir = Path(
        "E:/220204_BQ_dataset/Stage_info")
    orig_features = pyr_calc_all_features(_dir).transpose()
    # just get the scans and tumour labels.

    logging.info("Writing Original Features")

    # split metadata into separate columns
    _metadata = orig_features.index.str.split('_', expand=True).to_frame(index=False,
                                                                         name=['Date', 'Exp', 'Contour_Method',
                                                                               'Tumour_Model', 'Position', 'Age',
                                                                               'Cage_No.', 'Animal_No.'])

    orig_features.reset_index(inplace=True, drop=True)
    orig_features = pd.concat([_metadata, orig_features], axis=1)

    orig_features.index.rename('scanID', inplace=True)

    orig_features.to_csv(str(_dir / "orig_features.csv"))

    # get the images and masks
    logging.info("Getting values from inside the stage")
    scans_imgs, scan_names, masks = get_images_from_masks(_dir)

    logging.info("Counting Voxels")

    int_norm = normalise.NonRegMaskNormalise()

    # normalise the images, use first scan as the reference ???

    logging.info(f"Using {scan_names[0]} as the reference")
    int_norm.add_reference(scans_imgs[0], masks[0])

    logging.info("Normalising")
    int_norm.normalise(scans_imgs, masks)

    logging.info('writing normalised files')

    for i, vol in enumerate(scans_imgs):
        normed_dir_ = str(_dir / 'normed' / scan_names[i].split("/")[-1]) + ".nrrd"
        nrrd.write(normed_dir_, vol)

    logging.info("Recalculating Features")
    normed_features = pyr_calc_all_features(_dir, normed=True).transpose()

    _metadata = normed_features.index.str.split('_', expand=True).to_frame(index=False,
                                                                           name=['Date', 'Exp', 'Contour_Method',
                                                                                 'Tumour_Model', 'Position', 'Age',
                                                                                 'Cage_No.', 'Animal_No.'])
    normed_features.reset_index(inplace=True, drop=True)
    normed_features = pd.concat([_metadata, normed_features], axis=1)

    normed_features.index.rename('scanID', inplace=True)

    normed_features.to_csv(str(_dir / "normed_features.csv"))

    normed_comparison = pd.concat([_metadata,
                                   normed_features.compare(orig_features, keep_shape=False).rename(
                                       columns=lambda s: s.replace(["other", "self"],["original", "normed"]),
                                   inplace=True)],
                                  axis=1)

    normed_comparison.to_csv(_dir / "norm_comparision.csv")
    logging.info("DONE")


if __name__ == '__main__':
    main()
