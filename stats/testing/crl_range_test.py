#!/usr/bin/env python

from os.path import join, basename
import os
import sys
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common
import logging
import csv
import yaml
import shutil
import subprocess as sub
from lama_stats import LamaStats
import numpy as np
import SimpleITK as sitk


crl_string = """20140122_SCN4A_18.1_e_wt_rec_scaled_3.1241_pixel_14.nrrd	0.86325675941
20140128_SMOC1_18.2_c_wt_rec_scaled_3.1241_pixel_14.nrrd	0.889313191768
20150409_JAG2_E14.5_13.3h_WT_XX_rec_scaled_4.6823_pixel_14.0001.nrrd	0.892344286105
20141031_EPAS1_E14.5_16.1b_WT_XY_REC_scaled_4.6878_pixel_14.nrrd	0.897932479897
20160603_KCNJ10_E18.5_3.1g_WT_XY_REC_scaled_4.7297_pixel_13.9999.nrrd	0.898210579664
20140212_OTUB1_16.1e_wt_rec_scaled_3.1241_pixel_14.nrrd	0.903502236016
20141031_EPAS1_E14.5_16.1d_WT_XX_REC_scaled_4.6878_pixel_14.nrrd	0.912799789586
20140430_KLF7_E14.5_16.5f_WT_XY_rec_scaled_2.8248_pixel_14.nrrd	0.925373853996
20140515_KLHDC2_E14.5_21.1h_WT_XX_REC_scaled_2.8248_pixel_14.nrrd	0.945459124246
20140121_RIC8B_15.5_h_wt_rec_scaled_3.1241_pixel_14.nrrd	0.947180302834
20131206_MLLT3_15.3_d_WT_rec_scaled_5.1546_pixel_14.nrrd	0.956627654488
20140226_FGF10_E14.5_13.2_d_WT_rec_scaled_3.1241_pixel_14.nrrd	0.95853434824
20141021_CTNNBIP1_E14.5_13.4e_WT_XX_REC_scaled_4.6878_pixel_14.nrrd	0.961260771144
20140910_F10_E14.5_17.2a_WT_XX_REC_scaled_2.8248_pixel_14.nrrd	0.96884456294
20150220_Gabarapl2_E14.5_18.1i_wt_xy_REC_scaled_4.6878_pixel_14.nrrd	0.977330417923
20141119_CASZ1_E14.5_22.2a_WT_XY_REC_scaled_4.6878_pixel_14.nrrd	0.980437603699
20150216_Otud6b_E14.5_19.4c_WT_xy_REC_scaled_4.6878_pixel_14.nrrd	0.992673574071
20140121RIC8B_15.4_b_wt_rec_scaled_3.1241_pixel_14.nrrd	0.995529786698
20141120_GABARAPL2_E14.5_16.3h_WT_XX_rec_scaled_4.6823_pixel_14.0001.nrrd	0.997513541249
20141021_CTNNBIP1_E14.5_13.4d_WT_XX_scaled_4.6878_pixel_14.nrrd	1.00118023936
20140430_KLF7_E14.5_16.5b_WT_XX_rec_scaled_2.8248_pixel_14.nrrd	1.00462210171
20140409_ATG3_E14.5_16.4i_WT_XX_REC_scaled_2.8248_pixel_14.nrrd	1.01026325031
20141120_GABARAPL2_E14.5_20.1i_WT_XY_REC_scaled_4.6878_pixel_14.nrrd	1.01544253203
20150526_NTRK1_E14.5_11.2b_WT_XX_rec_scaled_4.6878_pixel_14.0.nrrd	1.02058727452
20140212_OTUB1_16.3e_wt_rec_scaled_3.1241_pixel_14.nrrd	1.02118963618
20151127_SLC16A3_E14.5_10.4e_WT_XY_rec_scaled_4.6878_pixel_14.0.nrrd	1.02470128482
20140409_ATG3_E14.5_18.2h_WT_XY_REC_scaled_2.8248_pixel_14.nrrd	1.0269981065
20160406_ATP1A2_E14.5_2.3g_WT_ND_scaled_4.7297_pixel_13.9999.nrrd	1.03642691295
20160323_SYNJ1_E14.5_10.4f_WT_XX_rec_scaled_4.7297_pixel_13.9999.nrrd	1.04626830835
20141030_PBX3_E14.5_13.5a_WT_XY_rec_2_scaled_4.6823_pixel_14.0001.nrrd	1.04700321469
20150411_RAPSN_E14.5_16.1f_WT_XY_rec_scaled_4.6823_pixel_14.0001.nrrd	1.05088472091
20140815_NXN_E14.5_17.4g_WT_XY_rec_scaled_2.8248_pixel_14.nrrd	1.06778640916
20141029_ETHE1_E14.5_13.4f_WT_XX_REC_scaled_4.6878_pixel_14.nrrd	1.07595924778
20141127_PRKAB2_E14.5_17.3D_WT_XY_rec_scaled_4.6878_pixel_14.nrrd	1.08786121345
20151129_PCX_E14.5_12.1d_WT_XX_rec_scaled_4.6878_pixel_14.0.nrrd	1.11024046484
20160608_GNAO1_E14.5_15.4i_WT_XX_rec_scaled_4.7297_pixel_13.9999.nrrd	1.13463101685
20140107_MLLT3_16.3_e_wt_rec_scaled_5.1546_pixel_14.nrrd	1.17077917172
20160125_ATP1A3_E14.5_7.5b_WT_XX_rec_scaled_4.6823_pixel_14.0001.nrrd	1.17801658513
20150305_HHIPL1_E14.5_19.1h_WT_XY_REC_scaled_4.6878_pixel_14.0.nrrd	1.22501187679"""


def make_groups_subset_files(mut_ids, wt_ids, outpath_wt, outpath_mut):
    with open(outpath_mut, 'wb') as cw:
        for m in mut_ids:
            cw.write("{}\n".format(m.strip('.nrrd')))
    with open(outpath_wt, 'wb') as cw:
        for w in wt_ids:
            cw.write("{}\n".format(w.strip('.nrrd')))


crls = []
jacobian_dir = "/media/neil/2c0a0602-17c4-4e12-9d8c-7412668fb150/work/male_female_jacobians_192_12/192_to_12"
root_out_dir = "/media/neil/2c0a0602-17c4-4e12-9d8c-7412668fb150/work/male_female_jacobians_192_12/13_wts_crl_range"
log_file = join(root_out_dir, 'results')
mask_path = "/media/neil/2c0a0602-17c4-4e12-9d8c-7412668fb150/work/male_female_jacobians_192_12/otsu_binary_img.nrrd"
stats_config = "/media/neil/2c0a0602-17c4-4e12-9d8c-7412668fb150/work/male_female_jacobians_192_12/13_wts_crl_range/stats.yaml"


for line in crl_string.splitlines():
    id_, scale_factor = line.split()
    crls.append((id_, scale_factor))
rev_crls = list(reversed(crls))


N = 13
mutant_sizes = [3, 4, 5, 6, 7, 8]
target_size = 9.18

with open(log_file, 'w') as results:
    for ms in mutant_sizes:
        mut_size_directory = join(root_out_dir, str(ms))
        common.mkdir_force(mut_size_directory)
        for f in range(8):
            results.write("mutant n:{}, step {}\n".format(str(ms), str(f)))
            diffs = []
            hit_voxels = []
            small = set([x[0] for x in crls[f: f+N]])
            small_scales = [float(x[1]) * target_size for x in crls[f: f+N]]

            large = set([x[0] for x in rev_crls[f: f+ms]])
            large_scales = [float(x[1]) * target_size for x in rev_crls[f: f + ms]]
            # Stop when we have any of the same ids present in each
            if len(small.union(large)) < N + ms:
                break
            outdir = join(mut_size_directory, str(f))
            common.mkdir_force(outdir)
            config = join(outdir, 'stats.yaml')
            shutil.copy(stats_config, config)
            wt_subset_file = join(outdir, 'wt_subset.csv')
            mut_subset_file = join(outdir, 'mut_subset.csv')
            make_groups_subset_files(small, large, wt_subset_file, mut_subset_file)
            LamaStats(config)
            # calculate mean difference
            mean_small = np.mean(small_scales)
            mean_large = np.mean(large_scales)
            diff = mean_large - mean_small
            diffs.append(str(diff))
            # get the number of hit resutls
            tsats_path = join(outdir, "jacobians_192_to_12/LinearModelR/__jacobians_192_to_12_genotype_LinearModelR_genotype_FDR_0.5_stats_.nrrd")
            img = sitk.ReadImage(tsats_path)
            arr = sitk.GetArrayFromImage(img)
            hit_voxel = arr[arr != 0].size
            hit_voxels.append(str(hit_voxel))
        results.write(",".join(diffs))
        results.write('\n')
        results.write(",".join(hit_voxels))
        print(f, diffs, hit_voxels)





