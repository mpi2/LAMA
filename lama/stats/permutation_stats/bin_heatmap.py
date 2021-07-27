"""
bin heatmaps (i.e. jacobians and intensities) into in a 3 * 4 * 8 array for permutation testing
"""

from pathlib import Path
import nrrd
import SimpleITK as sitk
import toml
from scipy import ndimage
from lama import common
import pandas as pd
import numpy as np
import os
import argparse


def get_heatmaps(dir, s):
    heatmap_list = []
    spec_name_list = []
    map_paths = [spec_path for spec_path in common.get_file_paths(dir) if ('log_jac' in str(spec_path))]
    # enumerating for speed only
    print(map_paths)
    for i, heatmap_path in enumerate(map_paths):
        heatmap, map_h = nrrd.read(heatmap_path)
        # only get heatmap vals inside of the mask

        spec_name_list.append(os.path.splitext(heatmap_path.name)[0].replace("log_jac_",""))
        heatmap = heatmap[s[0].start:s[0].stop,
                  s[1].start:s[1].stop,
                  s[2].start:s[2].stop]
        heatmap_list.append(heatmap)
    return heatmap_list, spec_name_list


def make_blocks_vectorized(x, d):
    p, m, n = x.shape

    # extend each axis of the array until it is divisible by x
    while p % d != 0:
        p += 1
    while m % d != 0:
        m += 1
    while n % d != 0:
        n += 1

    # reshape the array
    ref_shape = np.zeros((p, m, n,))

    ref_shape[:x.shape[0], :x.shape[1], :x.shape[2]] = x

    x = ref_shape

    return x.reshape(-1, m // d, d, n // d, d).transpose(1, 3, 0, 2, 4).reshape(-1, d, d, d)

def main():
    wt_dir = Path(
        "Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210423_g_by_e_stand_out/210415_g_by_e_anal/g_by_e_data/baseline/output/baseline")
    mut_dir = Path(
        "Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210423_g_by_e_stand_out/210415_g_by_e_anal/g_by_e_data/mutants/output/mutants")
    mask, mask_h = nrrd.read(
        "Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210423_g_by_e_stand_out/210415_g_by_e_anal/target/stats_mask.nrrd")
    ref, ref_h = nrrd.read(
        "Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210423_g_by_e_stand_out/210415_g_by_e_anal/target/210224_pop_avg_deformable_8.nrrd")

    # get bounding box of mask for cropping to mask
    # used multiple times
    s = ndimage.find_objects(mask)[0]

    # get reference cubes
    ref = ref[s[0].start:s[0].stop,
          s[1].start:s[1].stop,
          s[2].start:s[2].stop]  # crop to mask

    # ref_binned = make_blocks_vectorized(ref, 40)

    # write the reference cubes:
    # for i, cube in enumerate(ref_binned):
    #    nrrd.write("ref_cube_" + str(i) + ".nrrd", cube)

    # get heatmaps
    wt_htmps, wt_names = get_heatmaps(wt_dir, s)

    # turn heatmap into rubik's cube
    #wt_arrays = []

    #for htmp in wt_htmps:
    #    binned = make_blocks_vectorized(htmp, 40)
    #    # Summarise each bin by the non-zero mean. i.e. the faces/stickers of the
    #    # Rubik's cube
    #    face_val = [np.mean(cube[cube != 0]) for cube in binned]

    #    wt_arrays.append(face_val)

    #write to csv
    wt_df = pd.DataFrame(wt_htmps, index=wt_names)
    wt_df.to_csv("test_wt.csv")

    mut_htmps, mut_names = get_heatmaps(mut_dir, s)

    # turn heatmap into rubik's cube
    #mut_arrays = []

    #for htmp in mut_htmps:
    #    binned = make_blocks_vectorized(htmp, 40)
        # Summarise each bin by the non-zero mean. i.e. the faces/stickers of the
        # Rubik's cube
    #    face_val = [np.mean(cube[cube != 0]) for cube in binned]

    #    mut_arrays.append(face_val)

    #write to csv
    wt_df = pd.DataFrame(mut_htmps, index=mut_names)
    wt_df.to_csv("test_mut.csv")

if __name__ == '__main__':
    main()
