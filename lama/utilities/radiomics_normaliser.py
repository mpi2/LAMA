"""Normalises the radiomics scans by the average intensity of a mask"""
from typing import Union

from pandas import Series, DataFrame

from lama.img_processing import normalise
from lama import common
from lama.utilities import lama_img_info
from logzero import logger as logging
from lama import common
import os
import nrrd
from pathlib import Path

import numpy as np
import SimpleITK as sitk

from radiomics import featureextractor

import pandas as pd

from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas


# each scan in Ben's dataset will need its own mask
def get_images_from_masks(dir):
    img_list = []
    spec_name_list = []
    mask_list = []
    scan_paths = [spec_path for spec_path in common.get_file_paths(dir) if ('imgs' in str(spec_path))]
    mask_paths = [mask_path for mask_path in common.get_file_paths(dir) if ('labels' in str(mask_path))]

    # enumerate for indexing masks
    for i, img_path in enumerate(scan_paths):
        loader = common.LoadImage(img_path)
        img = loader.img
        # get_arr
        # cant use sitk.readimage due to header?
        # img, img_head = nrrd.read(img_path)
        # img = sitk.GetImageFromArray(img)

        m_loader = common.LoadImage(mask_paths[i])
        mask = m_loader.img

        # Only get values inside of the mask
        logging.info(f"Obtaining values from {img_path}")

        # s = ndimage.find_objects(mask)[0]

        # get the arrays
        # img_a = sitk.GetArrayFromImage(img)

        # mask_a = sitk.GetArrayFromImage(mask)

        #
        # img_a[(mask_a != 1) | (img_a < 0)] = 0

        # img_pro = sitk.GetImageFromArray(img_a)
        # mask_pro = sitk.GetImageFromArray(mask_a)
        # img_pro.CopyInformation(img)
        # mask_pro.CopyInformation(mask)

        spec_name_list.append(os.path.splitext(img_path.name)[0])
        # print(spec_name_list)
        img_list.append(img)
        mask_list.append(mask)
    return img_list, spec_name_list, mask_list


def pyr_calc_all_features(dir, normed: bool = False, images: list = None, file_names: list = None):
    # get either the normalised or original images
    scan_paths = images if normed \
        else [spec_path for spec_path in common.get_file_paths(dir) if ('imgs' in str(spec_path))]

    tumour_paths = [spec_path for spec_path in common.get_file_paths(dir) if ('tumour_respaced' in str(spec_path))]

    # debugging - Thanks Neil
    scan_paths.sort()
    tumour_paths.sort()

    # Get the first order measurements
    full_orders = []

    for i, img_path in enumerate(scan_paths):
        #logging.info(f"Calculating for {os.path.splitext(os.path.basename(img_path))[0]}")
        if normed: #files exist
            img = img_path
        else:
            loader = common.LoadImage(img_path)
            img = loader.img
            

        m_loader = common.LoadImage(tumour_paths[i])
        mask = m_loader.img

        # get all features and append to list
        extractor = featureextractor.RadiomicsFeatureExtractor()
        result = extractor.execute(img, mask)

        if file_names is not None:
            first_orders = pd.DataFrame.from_dict(result, orient='index',
                                                 columns=[os.path.splitext(os.path.basename(file_names[i]))[0]])
        else:
            first_orders= pd.DataFrame.from_dict(result, orient='index',
                                              columns=[os.path.splitext(os.path.basename(img_path))[0]])
        full_orders.append(first_orders)

    # fixing data format
    features = pd.concat(full_orders, axis=1).transpose()
    _metadata = features.index.str.split('_', expand=True).to_frame(index=False,
                                                                    name=['Date', 'Exp', 'Contour_Method',
                                                                          'Tumour_Model', 'Position', 'Age',
                                                                          'Cage_No.', 'Animal_No.'])
    features.reset_index(inplace=True, drop=True)
    features = pd.concat([_metadata, features], axis=1)

    features.index.rename('scanID', inplace=True)

    return features


def pyr_normaliser(_dir, _normaliser, scans_imgs, masks, fold: bool = False):
    # create a copy so orginal files aren't overwritten
    scans_imgs = scans_imgs.copy()

    # Do the normalisation
    if isinstance(_normaliser, normalise.NonRegMaskNormalise):
        _normaliser.add_reference(scans_imgs[0], masks[0])
        _normaliser.normalise(scans_imgs, masks, fold=fold, temp_dir=_dir)
    elif isinstance(_normaliser, normalise.IntensityHistogramMatch):
        _normaliser.normalise(scans_imgs, scans_imgs[0])

    return scans_imgs


def main():
    logging.info("Calculating Original First Order Features")
    _dir = Path("E:/220204_BQ_dataset/220307_BQ_norm")

    orig_features = pyr_calc_all_features(_dir)
    orig_features.to_csv(str(_dir / "orig_features.csv"))

    # get the images and masks
    logging.info("Getting values from inside the stage")
    scans_imgs, scan_names, masks = get_images_from_masks(_dir)

    # Copy scans and maks
    # scans_imgs_fold = scans_imgs.copy()
    # masks_fold = masks.copy()

    logging.info("Normalising to mean of the stage (subtraction)")
    sub_int_normed = pyr_normaliser(_dir, normalise.NonRegMaskNormalise(), scans_imgs, masks,)
    logging.info("Recalculating Features")
    sub_normed_features = pyr_calc_all_features(_dir, normed=True, images=sub_int_normed, file_names=scan_names)
    sub_normed_features.to_csv(str(_dir / "fold_normed_features.csv"))

    logging.info("Normalising to mean of the stage (fold)")
    fold_int_normed = pyr_normaliser(_dir, normalise.NonRegMaskNormalise(), scans_imgs, masks, fold=True, file_names=scan_names)
    logging.info("Recalculating Features")
    fold_normed_features = pyr_calc_all_features(_dir, normed=True, images=fold_int_normed)
    fold_normed_features.to_csv(str(_dir / "fold_normed_features.csv"))

    logging.info("Maskless Histogram Intensity Matching")
    histo_normed = pyr_normaliser(_dir, normalise.IntensityHistogramMatch(), scans_imgs, masks, file_names=scan_names)
    logging.info("Recalculating Features")
    histo_normed_features = pyr_calc_all_features(_dir, normed=True, images=histo_normed)
    histo_normed_features.to_csv(str(_dir / "fold_normed_features.csv"))

    all_features = pd.concat([orig_features, sub_normed_features, fold_normed_features, histo_normed_features],
                             keys=["Raw", "Subtraction", "Fold", "Histogram"])

    all_features.to_csv(str(_dir / "all_features.csv"))

    logging.info("Plotting Results")


    plot = all_features.plot(kind='box', subplots=True)

    fig = plot.get_figure()

    fig.savefig(str(_dir / "feature_comparision.png"))

    # normed_comparison = pd.concat([_metadata,
    #                               normed_features.compare(orig_features, keep_shape=False)],
    #                              axis=0)

    # normed_comparison.columns = [''.join(col) for col in normed_comparison.columns]

    # TODO: make this less clunky
    # normed_comparison.rename(
    #    columns=lambda s: s.replace("other", "--raw"),
    #    inplace=True)

    # normed_comparison.rename(
    #    columns=lambda s: s.replace("self", "--norm_sub"),
    #    inplace=True)

    # normed_comparison.to_csv(_dir / "fold_norm_comparision.csv")
    logging.info("DONE")




if __name__ == '__main__':
    main()
