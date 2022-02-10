from lama.img_processing import normalise
from logzero import logger as logging
from lama import common
import os
import nrrd
from pathlib import Path
from scipy import ndimage

from lama.utilities.radiomics_normaliser import get_images_from_masks

def main():
    _dir = Path(
        "E:/220204_BQ_dataset/stage_info/")

    # get the images and masks
    out_dir = _dir / 'out_dir'

    scans_imgs, scan_names = get_images_from_masks(_dir)



    # int_norm = normalise.IntensityMaskNormalise()

    # normalise the images
    # int_norm.add_reference(scans_imgs)

    #int_norm.normalise(mut_imgs)

    #int_norm.normalise(scans_imgs)

    #logging.info('writing normalised files')


    logging.info(f" saving files to {out_dir}")
    for i, vol in enumerate(scans_imgs):
        file_name = out_dir / scan_names[i]
        print(file_name)
        nrrd.write(file_name +".nrrd", vol)


if __name__ == '__main__':
    main()