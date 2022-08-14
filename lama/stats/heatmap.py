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


def heatmap(data: pd.DataFrame, title, use_sns=False):
    print("heatmap data: ", data)
    fig, ax = plt.subplots(figsize=[14, 15])
    # use_sns = False
    if use_sns:
        # sns.palplot(sns.color_palette("coolwarm"))
        if data.isnull().values.all():
            return
        try:
            sns.heatmap(data, center = 1.00, cmap=sns.color_palette("coolwarm", 100), ax=ax, cbar_kws={'label': 'mean volume ratio'},
                        square=True)
        except ValueError:
            ...

        ax.figure.axes[-1].yaxis.label.set_size(22)
        cbar = ax.collections[0].colorbar
        cbar.ax.tick_params(labelsize=20)
    else:
        ax.imshow(data)

    xlabels = data.columns
    ylabels = [x.replace('_', ' ') for x in data.index]
    # We want to show all ticks...
    # ax.set_xticks(np.arange(len(xlabels)))
    # ax.set_yticks(np.arange(len(ylabels)))
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
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    return True


def clustermap(data: pd.DataFrame, title, use_sns=False):
    print("heatmap data: ", data)
    print(np.std(data))
    fig, ax = plt.subplots(figsize=[14, 15])
    # use_sns = False
    if use_sns:
        # sns.palplot(sns.color_palette("coolwarm"))
        if data.isnull().values.all():
            return
        try:

            sns.clustermap(data,
                           z_score=0,
                           metric="euclidean",
                           cmap=sns.diverging_palette(250, 15, l=70, s=400,sep=40, n=512, center="light", as_cmap=True),
                           cbar_kws={'label': 'mean volume ratio'},
                           square=True)
        except ValueError as e:
            print(e)

            ...

        # ax.figure.axes[-1].yaxis.label.set_size(22)
        # cbar = ax.collections[0].colorbar
        # cbar.ax.tick_params(labelsize=20)
    else:
        ax.imshow(data)

    # xlabels = data.columns
    # ylabels = [x.replace('_', ' ') for x in data.index]
    # We want to show all ticks...
    # ax.set_xticks(np.arange(len(xlabels)))
    # ax.set_yticks(np.arange(len(ylabels)))
    # ... and label them with the respective list entries
    # ax.set_xticklabels(xlabels, rotation = 90, ha="center")
    # ax.set_yticklabels(ylabels)

    # Rotate the tick labels and set their alignment.
    # plt.setp(ax.get_xticklabels(), rotation=90, ha='center',
    #          rotation_mode="default")
    # ax.xaxis.set_tick_params(horizontalalignment='left')
    # plt.xticks(np.arange(len(gamma_range)) + 0.5, gamma_range, rotation=45, )
    # plt.xticks(rotation=90)
    # ax.tick_params(axis="x", direction="right", pad=-22)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)

    return True
