from logzero import logger as logging
from lama import common
import os
import nrrd
from pathlib import Path

import numpy as np
import SimpleITK as sitk

import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

pdf = PdfPages("E:/220204_BQ_Dataset/220508_BQ_norm/feature_comparision.pdf")


def multiple_plot(data):
    plt.figure()
    plt.clf()

    logging.info(data.columns[4])

    plt.xlabel('')
    # plt.plot(data)
    #data.groupby(['original_shape_VoxelVolume']).apply(print)

    #data.groupby(['Norm_Type']).plot.scatter(
    #    x='original_shape_VoxelVolume',
    #    y=2,
    #    hue='Norm_Type',
    #    #subplots=True,
    #    by=['Norm_Type'],
    #    #layout=(1, 4),
    #    sharex=False,
    #    rot=90,
    #    fontsize=5
    #)

    fig, ax = plt.subplots(figsize=(8, 6))

    #sns.scatterplot(x='original_shape_VoxelVolume',y=str(data.columns[3]),
    #            hue='Exp',
    #            data=data)


    data.groupby(['Norm_Type']).boxplot(
        by=['Exp', 'Age', 'Tumour_Model'],
        layout=(1, 4),
        sharex=True,
        rot=45,
        fontsize=5
    )

    plt.suptitle(str(data.columns[4]),fontsize='medium')

    plt.xlabel('')
    # plt.show()
    pdf.savefig()
    # plt.close()


def main():
    all_features = pd.read_csv("E:/220204_BQ_Dataset/220508_BQ_norm/all_features.csv")

    # all_features = all_features.pivot(index="Norm_type", columns='scanID')

    # print(all_features.iloc[32:np.shape(all_features)[1]].iloc[:, 0])

    data = all_features.iloc[:, 32:np.shape(all_features)[1]]

    logging.info("Plotting Results")

    # fig, ax = plt.subplots(len(all_features.columns), 1)

    for i, col in enumerate(data):
        plot_data = pd.concat([all_features['Norm_Type'],
                              #all_features['ScanID'],
                               #all_features['original_shape_VoxelVolume'],
                               all_features['Exp'],
                               all_features['Tumour_Model'],
                               all_features['Age'],
                               data.iloc[:, i]], axis=1)
        if data.iloc[:, i].name != 'original_shape_VoxelVolume':
            multiple_plot(plot_data)

    # all_features.plot(x='Norm_type',
    #                  kind='box',
    #                  subplots=True,
    #                  layout=(len(all_features.columns), 1),
    #                  legend=True,
    #                  figsize=(29.7, 21.1))
    # plt.legend(loc='best')

    # plt.savefig(str("E:/220204_BQ_Dataset/220331_BQ_norm/feature_comparision.pdf"))
    plt.close()
    pdf.close()


if __name__ == '__main__':
    main()
