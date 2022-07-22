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

#pdf = PdfPages("E:/220617_BQ_norm_body_full/feature_comparision.pdf")
#pdf = PdfPages("E:/Bl6_data/220524_test_radiomics/feature_comparision.pdf")

def multiple_plot(data):
    plt.figure()
    plt.clf()

    logging.info(data.columns[4])

    plt.xlabel('')
    #plt.plot(data)
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


    #sns.scatterplot(x='original_shape_VoxelVolume',y=str(data.columns[4]),
    #            hue='Genotype',
    #            style='Norm_Type',
    #            data=data)

    fig, ax = plt.subplots(figsize=(8, 6))



    data.groupby(['Tumour_Model']).boxplot(
        by=['Exp', 'Age', 'Tumour_Model'],
        #by=['Genotype','Background'],
        layout=(1, 4),
        sharex=True,
        rot=45,
        fontsize=5
    )

    plt.suptitle(str(data.columns[4]),fontsize='medium')

    plt.xlabel('')


    #plt.show()
    #pdf.savefig()
    plt.close()


def main():
    all_features = pd.read_csv("E:/220719_stage_norm_with_filts/sub_normed_features.csv", index_col=False)

    print(all_features)

    #all_features = pd.read_csv("E:/220508_BQ_norm/all_features.csv")

    #all_features = pd.read_csv("E:/220204_BQ_Dataset/220530_BQ_norm/histo_normed_features.csv")

    #all_features = pd.read_csv("E:/Bl6_data/220524_test_radiomics/lat_vent_features.csv")

    #all_features = all_features.pivot(index="Norm_Type", columns='scanID')


    #data = all_features.iloc[:, 32:np.shape(all_features)[1]]



    data = all_features.iloc[:, 32:np.shape(all_features)[1]]

    print("data", data)

    logging.info("Plotting Results")

    #fig, ax = plt.subplots(len(all_features.columns), 1)



    for i, col in enumerate(data):

        plot_data = pd.concat([#all_features['Norm_Type'],
                                all_features['scanID'],
                                all_features['Exp'],
                                all_features['Tumour_Model'],
                                all_features['Age'],
                                data.iloc[:, i]], axis=1)
         #if data.iloc[:, i].name != 'original_shape_VoxelVolume':
        #multiple_plot(plot_data)

        #plot_data = pd.concat([all_features['Norm_Type'],
                               #all_features['original_shape_VoxelVolume'],
        #                       all_features['Background'],
        #                       all_features['Genotype'],
        #                       data.iloc[:, i]], axis=1)
        #if data.iloc[:, i].name != 'original_shape_VoxelVolume':
        #multiple_plot(plot_data)

    # all_features.plot(x='Norm_type',
    #                  kind='box',
    #                  subplots=True,
    #                  layout=(len(all_features.columns), 1),
    #                  legend=True,
    #                  figsize=(29.7, 21.1))
    # plt.legend(loc='best')


    plt.close()
    #pdf.close()



    data = data[(all_features['Exp'] == 'MPTLVo7') &
                (all_features['Age'] == 'D14')]

    print(data)







    #data = data.set_index(all_features[all_features['Norm_Type'] == 'Subtraction']['Embryo'])
    data = data.set_index(all_features[(all_features['Exp'] == 'MPTLVo7') &
                                      (all_features['Age'] == 'D14')]['Tumour_Model'])
    #data = data.set_index(all_features['Tumour_Model'])
    print(data.index)
    print(data)



    #data.columns = data.columns.str.replace("original_",'')
    data = data.transpose()

    data = data.apply(lambda x: (x-x.mean())/x.std(), axis=1)

    data = data.apply(lambda x: (x - x.mean()) / x.std(), axis=0)


    #Drop na-cols
    data.dropna(axis='rows', inplace=True)



    fig, ax = plt.subplots(figsize=[56, 60])
    print(data)
    sns.clustermap(data,
                   figsize=[21, 21],
                   dendrogram_ratio=0.1,
                   metric="correlation",
                   #cmap=sns.diverging_palette(250, 15, l=70, s=400, sep=40, n=512, center="light", as_cmap=True),
                   #cbar_kws={'Genotype': 'Background'},
                   square=True,
                   xticklabels=True,
                   yticklabels=True)
    plt.tight_layout()

    plt.savefig("E:/220719_stage_norm_with_filts/radiomics_clustermap.png")
    plt.close()


if __name__ == '__main__':
    main()
