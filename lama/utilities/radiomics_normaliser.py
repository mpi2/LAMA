"""Normalises the radiomics scans by the average intensity of a mask"""
from typing import Union



from lama.img_processing import normalise

from logzero import logger as logging
from lama import common
import os
import nrrd
from pathlib import Path

import numpy as np
import SimpleITK as sitk

from radiomics import featureextractor, imageoperations
from scipy import ndimage
import pandas as pd

import raster_geometry as rg

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

        s = ndimage.find_objects(mask)[0]

        #get the arrays
        img_a = sitk.GetArrayFromImage(img)

        mask_a = sitk.GetArrayFromImage(mask)

        #
        img_a[(mask_a != 1) | (img_a < 0)] = 0

        # img_pro = sitk.GetImageFromArray(img_a)
        # mask_pro = sitk.GetImageFromArray(mask_a)
        # img_pro.CopyInformation(img)
        # mask_pro.CopyInformation(mask)

        spec_name_list.append(os.path.splitext(img_path.name)[0])
        # print(spec_name_list)
        img_list.append(img)
        mask_list.append(mask)
    return img_list, spec_name_list, mask_list

def spherify(dir):

    scan_paths = [spec_path for spec_path in common.get_file_paths(dir) if ('imgs' in str(spec_path))]
    tumour_paths =[spec_path for spec_path in common.get_file_paths(dir) if ('tumour_respaced' in str(spec_path))]


    # debugging - Thanks Neil
    scan_paths.sort()
    tumour_paths.sort()

    img_list = []
    spec_name_list = []
    mask_list = []

    for i, img_path in enumerate(scan_paths):

        # logging.info(f"Calculating for {os.path.splitext(os.path.basename(img_path))[0]}")

        logging.info(img_path)
        logging.info(tumour_paths[i])
        loader = common.LoadImage(img_path)
        img = loader.img

        m_loader = common.LoadImage(tumour_paths[i])
        mask = m_loader.img

        m_array = sitk.GetArrayFromImage(mask)

        s = ndimage.find_objects(m_array)[-1]



        midpoint = [np.round(np.mean([s[0].start, s[0].stop]))/512,
                    np.round((np.mean([s[1].start, s[1].stop]))) / 512,
                    np.round(482-(np.mean([s[2].start, s[2].stop]))) / 512]
        #print("Original Midpoint", [i*512 for i in midpoint])

        #print("Modified midpoint", midpoint)

        arr = rg.sphere(512, 10, midpoint, smoothing=True).astype(np.int_)


        ball = sitk.GetImageFromArray(arr)

        ball.CopyInformation(mask)

        sphere = "sphere"
        sitk.WriteImage(ball,
                        str(Path(os.path.dirname(img_path)).parent.absolute() / sphere/ os.path.basename(img_path)))

        spec_name_list.append(os.path.splitext(img_path.name)[0])

        # print(spec_name_list)
        #img_list.append(img)
        #mask_list.append(ball)

    return img_list, mask_list, spec_name_list

def pyr_calc_all_features(dir, normed: bool = False, images: list = None, file_names: list = None, spheres: list = None):
    # get either the normalised or original images
    scan_paths = images if normed else [spec_path for spec_path in common.get_file_paths(dir) if ('imgs' in str(spec_path))]

    tumour_paths = spheres if spheres\
        else [spec_path for spec_path in common.get_file_paths(dir) if ('tumour_respaced' in str(spec_path))]

    # debugging - Thanks Neil
    if not normed:
        scan_paths.sort()
    if not spheres:
        tumour_paths.sort()

    # Get the first order measurements
    full_orders = []

    for i, img_path in enumerate(scan_paths):

        #logging.info(f"Calculating for {os.path.splitext(os.path.basename(img_path))[0]}")
        if normed: #files exist
            img = img_path

        else:
            logging.info(img_path)
            logging.info(tumour_paths[i])
            loader = common.LoadImage(img_path)
            img = loader.img
        if spheres:
            mask = tumour_paths[i]
        else:
            m_loader = common.LoadImage(tumour_paths[i])
            mask = m_loader.img

        # apply pyradiomic filters
        #print(globals().items)
        #img_filt_ops = [x for x, y in globals().items if (x.startswith('pyradiomics.imageoperations.get') and x.endswith('Image'))]

        #print(img_filt_ops)

        # get all features and append to list
        extractor = featureextractor.RadiomicsFeatureExtractor(normalize=True)

        extractor.enableAllImageTypes()
        extractor
        extractor.enableAllFeatures()
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
    _metadata.reset_index(inplace=True, drop=True)
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
        print(type(scans_imgs[0]))
        _normaliser.normalise(scans_imgs, scans_imgs[0])

    return scans_imgs


def main():
    #import argparse
    #parser = argparse.ArgumentParser("Run various intensity normalisation methods")
    #parser.add_argument('-i', dest='indirs', help='dir with vols, tumour masks and label masks',
    #                    required=True)
    _dir = Path("E:/220204_BQ_dataset/220521_BQ_norm")
    #args = parser.parse_args()
    logging.info("Create Spheres from midpoint of tumour")

    images, spheres, scan_names = spherify(_dir)



    logging.info("Calculating Original Features")
    #_dir = Path(args.indirs)

    #orig_features = pyr_calc_all_features(_dir, spheres=spheres)
    #orig_features.to_csv(str(_dir / "orig_features.csv"))

    # get the images and masks
    logging.info("Getting values from inside the stage")
    scans_imgs, scan_names, masks = get_images_from_masks(_dir)

    scan_names.sort()
    logging.info("Normalising to mean of the stage (subtraction)")
    #sub_int_normed = pyr_normaliser(_dir, normalise.NonRegMaskNormalise(), scans_imgs, masks)



    #for i, vol in enumerate(sub_int_normed):
    #    file_name = scan_names[i] + '.nrrd'
    #    sitk.WriteImage(vol, str(_dir / file_name))
    #logging.info("Recalculating Features")
    #sub_normed_features = pyr_calc_all_features(_dir, normed=True, images=sub_int_normed, file_names=scan_names, spheres=spheres)
    #sub_normed_features.to_csv(str(_dir / "sub_normed_features.csv"))

    #logging.info("Normalising to mean of the stage (fold)")
    #fold_int_normed = pyr_normaliser(_dir, normalise.NonRegMaskNormalise(), scans_imgs, masks, fold=True)
    logging.info("Recalculating Features")
    #fold_normed_features = pyr_calc_all_features(_dir, normed=True, images=fold_int_normed, file_names=scan_names, spheres=spheres)
    #fold_normed_features.to_csv(str(_dir / "fold_normed_features.csv"))

    logging.info("Maskless Histogram Intensity Matching")
    histo_normed = pyr_normaliser(_dir, normalise.IntensityHistogramMatch(), scans_imgs, masks)
    logging.info("Recalculating Features")
    histo_normed_features = pyr_calc_all_features(_dir, normed=True, images=histo_normed, file_names=scan_names, spheres=spheres)
    histo_normed_features.to_csv(str(_dir / "fold_normed_features.csv"))

    #all_features = pd.concat([orig_features, sub_normed_features, fold_normed_features, histo_normed_features],
    #                         keys=["Raw", "Subtraction", "Fold", "Histogram"])

    #all_features.index.rename('Norm_Type', inplace=True)

    #all_features.to_csv(str(_dir / "all_features.csv"))

    logging.info("DONE")


if __name__ == '__main__':
    main()
