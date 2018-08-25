import sys
from os.path import join
import os
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common
import SimpleITK as sitk
import os
from os.path import join, basename
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import matplotlib
from skimage import measure

def cluster_images(paths, mask, outdir):
    imgs = []
    names = []
    mask_img = sitk.ReadImage(mask)
    mask_arr = sitk.GetArrayFromImage(mask_img)

    for i, vp in enumerate(paths):
        bn = basename(vp)
        names.append(bn)
        img = sitk.ReadImage(vp)
        arr = sitk.GetArrayFromImage(img)
        imgs.append(arr)

    # Now get pairwise dice coefficent for each image

    results = []
    for name, img1 in zip(names, imgs):
        metrics = []
        for img2 in imgs:
            dis = measure.compare_nrmse(img1, img2)

            metrics.append(dis)
        results.append(metrics)
    for i, name in enumerate(names):
        print(i, name)

    matplotlib.rcParams['lines.linewidth'] = 5.0
    Z = hierarchy.linkage(results, 'single')
    plt.figure()
    dn = hierarchy.dendrogram(Z)

    hierarchy.set_link_color_palette(['m', 'c', 'y', 'k'])
    hierarchy.set_link_color_palette(None)  # reset to default after use
    plt.show()



# in_dir1 = '/home/neil/sig/LAMA_results/E14.5/compare_cbx2_male_female_variances/male_female_and_cbx2/all_inputs/mutant_runs/ATP2a1/full_run_14um/stats/intensity/normalised_images/mutant'
# in_dir2 = '/media/neil/2c0a0602-17c4-4e12-9d8c-7412668fb150/work/test_jac_clustering/mut'
mask = '/home/neil/sig/LAMA_results/E14.5/compare_cbx2_male_female_variances/male_female_and_cbx2/all_inputs/full_run_14um/padded_target/mask_cleaned.nrrd'
out_dir = '/home/neil/temp'

indir = '/home/neil/sig/LAMA_results/E14.5/compare_cbx2_male_female_variances/male_female_and_cbx2/all_inputs/mutant_runs/ATP2a1/full_run_14um/stats/intensity_all_mutants/n1'
paths = [join(indir, x) for x in os.listdir(indir)]

jaccard(paths, mask, out_dir)
