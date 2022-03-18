from logzero import logger as logging
from lama import common
import os
import nrrd
from pathlib import Path

import numpy as np
import SimpleITK as sitk

import matplotlib.pyplot as plt

import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

pdf = PdfPages("E:/220204_BQ_Dataset/220307_BQ_norm/feature_comparision.pdf")


def multiple_plot(data):
    plt.figure()
    plt.clf()

    logging.info(data.columns[2])

    plt.xlabel('')
    #plt.plot(data)

    data.groupby(['Norm_type']).boxplot(
                                by=['Exp'],
                                layout=(1,4),
                                sharex=False,
                                rot=30
                                )

    plt.suptitle(str(data.columns[2]))

    plt.xlabel('')
    #plt.show()
    pdf.savefig()
    #plt.close()


def main():
    all_features = pd.read_csv("E:/220204_BQ_Dataset/220307_BQ_norm/all_features.csv")

    # all_features = all_features.pivot(index="Norm_type", columns='scanID')

    # print(all_features.iloc[32:np.shape(all_features)[1]].iloc[:, 0])

    data = all_features.iloc[:, 32:np.shape(all_features)[1]]



    logging.info("Plotting Results")

    # fig, ax = plt.subplots(len(all_features.columns), 1)



    for i, col in enumerate(data):

        plot_data = pd.concat([all_features['Norm_type'],
                               all_features['Exp'],
                               #all_features['Tumour_Model'],
                               #all_features['Age'],
                               np.log(data.iloc[:, i])], axis=1)
        multiple_plot(plot_data)



    # all_features.plot(x='Norm_type',
    #                  kind='box',
    #                  subplots=True,
    #                  layout=(len(all_features.columns), 1),
    #                  legend=True,
    #                  figsize=(29.7, 21.1))
    # plt.legend(loc='best')

    #plt.savefig(str("E:/220204_BQ_Dataset/220307_BQ_norm/feature_comparision.pdf"))
    plt.close()
    pdf.close()




if __name__ == '__main__':
    main()
