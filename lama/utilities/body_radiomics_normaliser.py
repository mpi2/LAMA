"""Normalises the radiomics scans by the average intensity of a mask"""
from typing import Union

from lama.img_processing import normalise

from logzero import logger as logging
from lama import common
import os
import nrrd
from pathlib import Path

import numpy as np
from scipy import ndimage
import SimpleITK as sitk

import radiomics

from radiomics import featureextractor

import pandas as pd


# each scan in Ben's dataset will need its own mask
def get_images_from_masks(dir):
    img_list = []
    spec_name_list = []
    mask_list = []
    scan_paths = [spec_path for spec_path in common.get_file_paths(dir) if ('imgs' in str(spec_path))]
    mask_paths = [mask_path for mask_path in common.get_file_paths(dir) if ('labels' in str(mask_path))]

    scan_paths.sort()
    mask_paths.sort()

    # enumerate for indexing masks
    for i, img_path in enumerate(scan_paths):

        # load image and stage
        logging.info(f"Obtaining values from {img_path}")
        loader = common.LoadImage(img_path)
        img = loader.img

        logging.info(f"Removing values from {mask_paths[i]}")
        m_loader = common.LoadImage(mask_paths[i])
        mask = m_loader.img

        # Only get values inside of the mask

        # get the arrays
        img_a = sitk.GetArrayFromImage(img)
        mask_a = sitk.GetArrayFromImage(mask)

        # remove the stage
        img_a[mask_a == 1] = np.min(img_a)

        img_pro = sitk.GetImageFromArray(img_a)
        img_pro.CopyInformation(img)

        logging.info(f"Performing Otsu on {scan_paths[i]}")

        Otsu = sitk.OtsuThresholdImageFilter()

        inv_mask = Otsu.Execute(img)
        o_mask = sitk.InvertIntensity(inv_mask, 1)

        o_mask = sitk.ConnectedComponent(o_mask != o_mask[0, 0, 0])

        # sitk.WriteImage(seg, os.path.join(output, name + "_all_connected.nrrd"))
        o_mask = sitk.RelabelComponent(o_mask)
        o_mask = o_mask == 1
        # sitk.WriteImage(seg, os.path.join(output, name + "_largest_connected.nrrd"))

        # lets see if dilate with a tight kernal fixes getting stupid dots everywhere.
        dilate = sitk.BinaryDilateImageFilter()
        dilate.SetKernelRadius([1, 1, 1])
        dilate.SetKernelType(sitk.sitkBall)
        o_mask = dilate.Execute(o_mask)
        o_mask.CopyInformation(img)

        # Operation is peformed using scipy so needs to be a numpy array
        npa = sitk.GetArrayFromImage(o_mask)

        logging.info("fill holes in first orientation")
        npa_hole_filled = fill_image(npa)

        logging.info("fill holes in second orientation")
        npa_hole_filled = fill_image(npa_hole_filled, roll=1)

        logging.info("fill holes in third orientation")
        npa_hole_filled = fill_image(npa_hole_filled, roll=0)

        # Need to turn the image back into its original orientation
        transposed = np.transpose(npa_hole_filled, axes=(0, 2, 1))

        # Turn np array to image
        filled = sitk.GetImageFromArray(transposed)
        filled.CopyInformation(o_mask)


        q_otsu = "otsu"
        sitk.WriteImage(filled, str(Path(os.path.dirname(img_path)).parent.absolute() / q_otsu / os.path.basename(img_path)))

        logging.info("Removing values outside of the mask")

        mask_a2 = sitk.GetArrayFromImage(filled)

        img_a2 = sitk.GetArrayFromImage(img_pro)

        img_a2[mask_a2 != 1] = np.min(img_a2)

        img_pro = sitk.GetImageFromArray(img_a2)
        img_pro.CopyInformation(img)

        spec_name_list.append(os.path.splitext(img_path.name)[0])
        # print(spec_name_list)
        img_list.append(img_pro)
        mask_list.append(filled)
    return img_list, spec_name_list, mask_list


def fill_image(npa, roll=0):
    """ Binary hole fill in any orientation

    Go slice by slice and fill any holes of a 3D image. Orientation is chosen by the roll parameter.

    :param array npa: Numpy array of the image to be processed
    :param int roll: Which axis to image is to be "rolled" in. If use none, 0 and 1. Should do all axis for a 3D image
                    [Default 0]
    :return array npa_hole_filled: Numpy array of image
    """

    npa_hole_filled = None

    # Change orientation
    if roll:
        loop = np.rollaxis(npa, roll)
    else:
        loop = npa

    # loop through each slice
    for slice_ in loop:
        slice_fill = ndimage.binary_fill_holes(slice_).astype(int)
        if npa_hole_filled is None:
            npa_hole_filled = slice_fill
        else:
            npa_hole_filled = np.dstack((npa_hole_filled, slice_fill))

    return npa_hole_filled


def pyr_calc_all_features(dir, normed: bool = False, images: list = None, file_names: list = None):
    # get either the normalised or original images
    scan_paths = images if normed \
        else [spec_path for spec_path in common.get_file_paths(dir) if ('imgs' in str(spec_path))]
    print(scan_paths)
    tumour_paths = [spec_path for spec_path in common.get_file_paths(dir) if ('tumour_respaced' in str(spec_path))]

    # debugging - Thanks Neil
    scan_paths.sort()
    tumour_paths.sort()

    # Get the first order measuremients
    full_orders = []

    for i, img_path in enumerate(scan_paths):

        # logging.info(f"Calculating for {os.path.splitext(os.path.basename(img_path))[0]}")
        if normed:  # files exist
            img = img_path
        else:
            logging.info(img_path)
            logging.info(tumour_paths[i])
            loader = common.LoadImage(img_path)
            img = loader.img

        m_loader = common.LoadImage(tumour_paths[i])
        mask = m_loader.img

        #so need to binarise the mask

        mask_arr = sitk.GetArrayFromImage(mask)

        mask_arr[mask_arr > 1] = 1

        b_mask = sitk.GetImageFromArray(mask_arr)
        b_mask.CopyInformation(mask)

        # get all features and append to list
        extractor = featureextractor.RadiomicsFeatureExtractor()


        result = extractor.execute(img, b_mask)

        if file_names is not None:
            first_orders = pd.DataFrame.from_dict(result, orient='index',
                                                  columns=[os.path.splitext(os.path.basename(file_names[i]))[0]])
        else:
            first_orders = pd.DataFrame.from_dict(result, orient='index',
                                                  columns=[os.path.splitext(os.path.basename(img_path))[0]])
        full_orders.append(first_orders)

    # fixing data format
    features = pd.concat(full_orders, axis=1).transpose()

    #_metadata = features.index #.str.split('_', expand=True).to_frame(index=False,
    #                                                                name=['Date', 'Strain', 'Colony',
    #                                                                      'Embryo', 'Genotype'])

                                                                    #name=['Date', 'Exp', 'Contour_Method',
                                                                    #      'Tumour_Model', 'Position', 'Age',
                                                                    #      'Cage_No.', 'Animal_No.'])
    #_metadata.reset_index(inplace=True, drop=True)
    #features.reset_index(inplace=True, drop=True)
    #features = pd.concat([_metadata, features], axis=1)

    #features.index.rename('scanID', inplace=True)

    return features


def pyr_normaliser(_dir, _normaliser, scans_imgs, masks, fold: bool = False, ref_vol_path: Path = None):
    # create a copy so orginal files aren't overwritten
    scans_imgs = scans_imgs.copy()

    # Do the normalisation
    if isinstance(_normaliser, normalise.NonRegMaskNormalise):
        _normaliser.add_reference(scans_imgs[0], masks[0])
        _normaliser.normalise(scans_imgs, masks, fold=fold, temp_dir=_dir)
    elif isinstance(_normaliser, normalise.IntensityHistogramMatch):
        if ref_vol_path:
            ref_vol = common.LoadImage(ref_vol_path).img
            _normaliser.normalise(scans_imgs, ref_vol)

        else:
            _normaliser.normalise(scans_imgs, scans_imgs[0])

    return scans_imgs


def main():
    # import argparse
    # parser = argparse.ArgumentParser("Run various intensity normalisation methods")
    # parser.add_argument('-i', dest='indirs', help='dir with vols, tumour masks and label masks',
    #                    required=True)

    # args = parser.parse_args()

    logging.info("Calculating Original Features")
    # _dir = Path(args.indirs)
    _dir = Path("E:/220204_BQ_dataset/220530_BQ_norm")

    #ref_path = Path("E:/Bl6_data/211014_g_by_back/target/210602_C3H_avg_n18.nrrd")

    #orig_features = pyr_calc_all_features(_dir)
    #orig_features.to_csv(str(_dir / "orig_features.csv"))

    # get the images and masks
    logging.info("Getting values from inside the stage")
    scans_imgs, scan_names, masks = get_images_from_masks(_dir)

    scan_names.sort()
#logging.info("Normalising to mean of the stage (subtraction)")
    
    sub_int_normed = pyr_normaliser(_dir, normalise.NonRegMaskNormalise(), scans_imgs, masks)

    # for i, vol in enumerate(sub_int_normed):
    #    file_name = scan_names[i] + '.nrrd'
    #    sitk.WriteImage(vol, str(_dir / file_name))
    # logging.info("Recalculating Features")
    sub_normed_features = pyr_calc_all_features(_dir, normed=True, images=sub_int_normed, file_names=scan_names)
    sub_normed_features.to_csv(str(_dir / "sub_normed_features.csv"))

    #logging.info("Normalising to mean of the stage (fold)")
    #fold_int_normed = pyr_normaliser(_dir, normalise.NonRegMaskNormalise(), scans_imgs, masks, fold=True)
    #logging.info("Recalculating Features")
    #fold_normed_features = pyr_calc_all_features(_dir, normed=True, images=fold_int_normed, file_names=scan_names)
    #fold_normed_features.to_csv(str(_dir / "fold_normed_features.csv"))

    #logging.info("Maskless Histogram Intensity Matching")
    #histo_normed = pyr_normaliser(_dir, normalise.IntensityHistogramMatch(), scans_imgs, masks, ref_vol_path=ref_path)
    #logging.info("Recalculating Features")
    #histo_normed_features = pyr_calc_all_features(_dir, normed=True, images=histo_normed, file_names=scan_names)
    #histo_normed_features.to_csv(str(_dir / "histo_normed_features.csv"))

    all_features = pd.concat([orig_features, sub_normed_features],
                             keys=["Raw", "Subtraction"])

    all_features.index.rename('Norm_Type', inplace=True)

    all_features.to_csv(str(_dir / "all_features.csv"))

    logging.info("DONE")


if __name__ == '__main__':
    main()
