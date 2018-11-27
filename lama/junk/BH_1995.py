#!/usr/bin/env python
import numpy as np
import SimpleITK as sitk

# This implementaion of benjamini Hochberg(1995) FDR correction gives almost exactly the same results as R's p.adjust

pvals_file ='/home/neil/sig/LAMA_results/E14.5/compare_cbx2_male_female_variances/male_female_and_cbx2/all_inputs/full_run_14um/stats_male_female/jacobians_male_female_192_to_12/LinearModelR/__jacobians_male_female_192_to_12_genotype_LinearModelR_pvals_genotype_stats_.nrrd'
tvals_file = '/home/neil/sig/LAMA_results/E14.5/compare_cbx2_male_female_variances/male_female_and_cbx2/all_inputs/full_run_14um/stats_male_female/jacobians_male_female_192_to_12/LinearModelR/__jacobians_male_female_192_to_12_genotype_LinearModelR_Tstats_genotype_stats_.nrrd'
mask_file = '/home/neil/sig/LAMA_results/E14.5/compare_cbx2_male_female_variances/male_female_and_cbx2/all_inputs/full_run_14um/padded_target/otsu_binary_img.nrrd'

Q = 0.05

full_pvals = sitk.GetArrayFromImage(sitk.ReadImage(pvals_file))
full_tvals = sitk.GetArrayFromImage(sitk.ReadImage(tvals_file))
mask = sitk.GetArrayFromImage(sitk.ReadImage(mask_file))

# Amsk the values
pvals = full_pvals[mask == 1]
tvals = full_tvals[mask == 1]

# Join the p and t valsue sinto a 2D array
p_t = np.dstack((pvals.flatten(), tvals.flatten()))[0]


# lowest_p = p_t[np.argmin(p_t[:, 0])]

# Sort the values bsaed on p-value
sorted_ = p_t[p_t[:, 0].argsort()]

# create a rank array
ranks = np.array(list(range(1, pvals.size + 1))).astype(np.float32)

# Gernate critical values
m = pvals.size
critical_values = (ranks/m) * Q

# Filter when above critiacla value
filtered_indices = np.where(sorted_[:, 0] < critical_values)
filtered_values = sorted_[filtered_indices]

# Get the t thresholds for positive and negative values
try:
    min_pos_t = filtered_values[:, 1][filtered_values[:, 1] > 0].min()
except ValueError:
    min_pos_t = None
try:
    max_neg_t = filtered_values[:, 1][filtered_values[:, 1] < 0].max()
except ValueError:
    max_neg_t = None

print(('min_pos_t', min_pos_t))
print(('max_neg_t', max_neg_t))




