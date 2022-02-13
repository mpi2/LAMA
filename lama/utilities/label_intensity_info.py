from lama.img_processing import normalise
from logzero import logger as logging
from lama import common
import os
import nrrd
from pathlib import Path
from scipy import ndimage
import numpy as np
import pandas as pd
from lama.utilities.radiomics_normaliser import get_images_from_masks
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt


def main():
    _dir = Path(
        "E:/220204_BQ_dataset/stage_info/out_dir")

    # get the images and masks
    # out_dir = _dir / 'out_dir'

    # scans_imgs, scan_names = get_images_from_masks(_dir)

    # int_norm = normalise.IntensityMaskNormalise()

    # normalise the images
    # int_norm.add_reference(scans_imgs)

    # int_norm.normalise(mut_imgs)

    # int_norm.normalise(scans_imgs)

    # logging.info('writing normalised files')

    # logging.info(f" saving masked files to {out_dir}")
    # for i, vol in enumerate(scans_imgs):
    #    file_name = out_dir / scan_names[i]

    #    nrrd.write(str(file_name) + ".nrrd", vol)

    vols = common.get_file_paths(_dir)

    all_vals = []
    for i, vol in enumerate(vols):
        img, img_h = nrrd.read(vol)
        all_vals.append(img[img != 0].flatten())

    info = pd.DataFrame([os.path.basename(vol) for vol in vols], columns=['name'])

    logging.info("basic stats")
    means = pd.DataFrame([np.mean(val) for val in all_vals])
    sds = pd.DataFrame([np.std(val) for val in all_vals])
    sems = pd.DataFrame([np.std(val) / np.sqrt(len(val)) for val in all_vals])

    bstats = pd.concat([info, means, sds, sems], axis=1)

    bstats.to_csv("E:/220204_BQ_dataset/stage_info/bstats.csv")
    d = pd.DataFrame(all_vals, index=info.values)

    df = pd.concat([info, d], axis=1)

    df = pd.melt(df, id_vars='name')

    logging.info("plotting boxplots")
    fig1, ax1 = plt.subplots()
    ax1.set_title('Basic Plot')
    ax1.boxplot(df)

    plt.show()

    logging.info("Stats")

    fit = smf.ols(formula='values ~ name', data=df, missing='drop').fit()
    p = fit.pvalues

    print(p)

    pd.DataFrame(p).to_csv("E:/220204_BQ_dataset/stage_info/pvals.csv")

    # all_vals = [vol[vol != 0].flatten() for vol in scans_imgs]

    # pd.DataFrame(all_vals).to_csv("E:/220204_BQ_dataset/stage_info/intens_vals.csv",
    #                              index=True)

    logging.info("Finished!")


if __name__ == '__main__':
    main()
