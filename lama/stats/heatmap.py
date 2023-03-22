"""
Seaborn 0.9.0 gives missing rows for some reason. So use matplotlib directly instead
"""

import matplotlib.pyplot as plt
import pandas as pd
from logzero import logger as logging
import numpy as np
import seaborn as sns
from typing import Union
import matplotlib
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from scipy.stats import zscore

def heatmap(data: pd.DataFrame, title, use_sns=False, rad_plot: bool = False):
    fig, ax = plt.subplots(figsize=[56, 60])
    # use_sns = False

    font_size = 14 if rad_plot else 22

    if use_sns:
        # sns.palplot(sns.color_palette("coolwarm"))
        if data.isnull().values.all():
            return
        try:
            sns.heatmap(data, center=1.00, cmap=sns.color_palette("coolwarm", 100), ax=ax,
                        cbar_kws={'label': 'mean volume ratio',
                                  'fraction': 0.05}, square=(not rad_plot))
        except ValueError:
            ...

        ax.figure.axes[-1].yaxis.label.set_size(font_size)
        cbar = ax.collections[0].colorbar
        cbar.ax.tick_params(labelsize=22)
        # adjust cbar fraction

    else:
        ax.imshow(data)

    xlabels = data.columns
    ylabels = [x.replace('_', ' ') for x in data.index]

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(xlabels)))

    ax.set_yticks(np.arange(len(ylabels)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(xlabels, rotation=90, ha="center")
    ax.set_yticklabels(ylabels)

    # Rotate the tick labels and set their alignment.
    # plt.setp(ax.get_xticklabels(), rotation=90, ha='center',
    #          rotation_mode="default")
    # ax.xaxis.set_tick_params(horizontalalignment='left')
    # plt.xticks(np.arange(len(gamma_range)) + 0.5, gamma_range, rotation=45, )
    # plt.xticks(rotation=90)
    # ax.tick_params(axis="x", direction="right", pad=-22)
    # Note for radiomics data make this stuff smaller
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)

    return True


def clustermap(data: pd.DataFrame, title, use_sns=False, rad_plot: bool = False):

    font_size = 10 if rad_plot else 22

    # use_sns = False
    if use_sns:
        # sns.palplot(sns.color_palette("coolwarm"))
        if data.isnull().values.all():
            return
        try:
            if rad_plot:
                cg = sns.clustermap(data,
                                metric="euclidean",
                                cmap=sns.diverging_palette(250, 15, l=70, s=400, sep=1, n=512, center="light",
                                                           as_cmap=True),
                                center=1,
                                cbar_kws={'label': 'mean volume ratio'}, square=True,
                                figsize=[30, len(data)*0.3])
            else:
                cg = sns.clustermap(data,
                                    z_score=0,
                                    metric="euclidean",
                                    cmap=sns.diverging_palette(250, 15, l=70, s=400, sep=40, n=512, center="light",
                                                               as_cmap=True),
                                    cbar_kws={'label': 'mean volume ratio'},
                                    figsize=[14, 15],
                                    square=True)

            ylabels = [x.replace('_', ' ') for x in data.index]

            cg.ax_heatmap.set_yticks(np.arange(len(ylabels)))
            cg.ax_heatmap.tick_params(axis='y', labelsize=font_size)
            cg.ax_heatmap.set_yticklabels(ylabels, fontsize=font_size, rotation=0)



        except ValueError as e:
            print(e)

            ...


    return True
