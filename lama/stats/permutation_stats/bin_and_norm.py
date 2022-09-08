from lama.img_processing import normalise
from logzero import logger as logging
from lama import common
import os
import nrrd
from pathlib import Path
from scipy import ndimage
import numpy as np
import SimpleITK as sitk

import pandas as pd

from lama.stats.permutation_stats import bin_heatmap
from lama.utilities import prep_for_man_valid as pfmv


def main():
    print("something")
    wt_dir = Path(
        "Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210519_int_anal/wt")
    mut_dir = Path(
        "Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210519_int_anal/non_wt")

    mask, mask_h = nrrd.read(
        "Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210423_g_by_e_stand_out/210415_g_by_e_anal/target/stats_mask.nrrd")

    pop_avg, pop_h = nrrd.read(
        "Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210423_g_by_e_stand_out/210415_g_by_e_anal/target/210224_pop_avg_deformable_8.nrrd")

    s = ndimage.find_objects(mask)[0]

    # get the images
    wt_imgs, wt_names = pfmv.get_images(wt_dir, s)

    mut_imgs, mut_names = pfmv.get_images(mut_dir, s)

    int_norm = normalise.IntensityMaskNormalise()

    # normalise the images
    int_norm.add_reference(wt_imgs)

    int_norm.normalise(mut_imgs)

    int_norm.normalise(wt_imgs)

    wt_arrays = []

    for img in wt_imgs:
        binned = bin_heatmap.make_blocks_vectorized(img, 40)
        # Summarise each bin by the non-zero mean. i.e. the faces/stickers of the
        # Rubik's cube
        face_val = [np.mean(cube[cube != 0]) for cube in binned]

        wt_arrays.append(face_val)

    # write to csv
    wt_df = pd.DataFrame(wt_arrays, index=wt_names)
    wt_df.to_csv("test_wt.csv")

    mut_arrays = []

    for img in mut_imgs:
        binned = bin_heatmap.make_blocks_vectorized(img, 40)
        # Summarise each bin by the non-zero mean. i.e. the faces/stickers of the
        # Rubik's cube
        face_val = [np.mean(cube[cube != 0]) for cube in binned]

        mut_arrays.append(face_val)

    # write to csv
    wt_df = pd.DataFrame(mut_arrays, index=mut_names)
    wt_df.to_csv("test_mut.csv")

if __name__ == '__main__':
    main()