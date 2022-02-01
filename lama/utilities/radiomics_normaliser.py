"""Normalises the radiomics scans by the average intensity of a mask"""
from lama.img_processing import normalise
from logzero import logger as logging
from lama import common
import os
import nrrd
from pathlib import Path
from scipy import ndimage
import numpy as np
import SimpleITK as sitk
from lama.utilities import man_clip_removal

# each scan in Ben's dataset will need its own mask
def get_images_and_masks(dir):
    img_list = []
    spec_name_list = []
    scan_paths = [spec_path for spec_path in common.get_file_paths(dir) if ('img' in str(spec_path))]
    mask_paths = [mask_path for mask_path in common.get_file_paths(dir) if ('better' in str(mask_path))]

    # enumerate for indexing masks
    for i, img_path in enumerate(scan_paths):
        print(mask_paths[i])
        mask, m_h = nrrd.read(mask_paths[i])
        print(mask)
        s = ndimage.find_objects(mask)[0]

        img, img_h = nrrd.read(img_path)
        img = img[s[0].start:s[0].stop,
              s[1].start:s[1].stop,
              s[2].start:s[2].stop]
        spec_name_list.append(os.path.splitext(img_path.name)[0])
        print(spec_name_list)
        img_list.append(img)
    return img_list, spec_name_list


def main():
    _dir = Path(
        "E:/BenQ_tumour_test")

    # get the images and masks

    scans_imgs, scan_names = get_images_and_masks(_dir)

    print(scan_names)

    int_norm = normalise.IntensityMaskNormalise()

    # normalise the images
    int_norm.add_reference(scans_imgs)

    #int_norm.normalise(mut_imgs)

    int_norm.normalise(scans_imgs)

    logging.info('writing normalised files')

    for i, vol in enumerate(scans_imgs):
        nrrd.write(scan_names[i] + "_normed.nrrd", vol)

    logging.info('writing normalised files')
if __name__ == '__main__':
    main()