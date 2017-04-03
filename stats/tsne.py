#!/usr/bin/env python

import SimpleITK as sitk
from os.path import join, basename
from sklearn.manifold import TSNE
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
from os.path import isdir
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common
from collections import OrderedDict


THRESH = 4.0
MAX_Z = 400  # If size extent is larger, then rescale to speed things up/reduce memory usage

# t-sne parameters
TSNE_PARAMETERS = {
    'n_components': 2,
    'init': 'pca',
    'random_state': 0,
    'n_iter': 60000,
    'perplexity': 5,
    'learning_rate': 6000,
    'method': 'barnes_hut',
    'metric': 'dice'
}


def cluster(indir, outpath):
    if isdir(indir):
        names = OrderedDict()
        paths = common.GetFilePaths(indir)
    else:
        paths = common.get_inputs_from_file_list(indir)

    imgs = []
    first_img = sitk.ReadImage(paths[0])
    if first_img.GetSize()[2] > MAX_Z:
        # determine scaling factor
        scale_factor = int(first_img.GetSize()[2] / MAX_Z)
        if scale_factor < 2:
            scale_factor = 2
    else:
        scale_factor = False

    for i, vp in enumerate(paths):
        print(vp)
        bn = basename(vp)
        names[i] = bn
        img = sitk.ReadImage(vp)
        if scale_factor:
            img = sitk.BinShrink(img, (scale_factor, scale_factor, scale_factor))
        arr = sitk.GetArrayFromImage(img)
        try:
            arr[np.where((arr > -THRESH) & (arr < THRESH))] = 0
        except ValueError:  # No values less than 4
            continue
        imgs.append(arr.ravel())

    tsne = TSNE(**TSNE_PARAMETERS)
    trans_data = tsne.fit_transform(imgs).T

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.scatter(trans_data[0], trans_data[1], cmap=plt.cm.rainbow)
    for i in range(trans_data[0].size):
        ax.annotate(names.keys()[i], xy=(trans_data[0][i], trans_data[1][i]))
    plt.title("t-SNE differences beteen mutants")

    plt.axis('tight')
    plt.savefig(outpath)

    return names


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("Create t-sne clustering plot of images")
    parser.add_argument('-i', '--indir', dest='indir', help='path to folder with images or a file with paths', required=True)
    parser.add_argument('-o', '--outpath', dest='outpath', help='path to file to save plot figure to', required=True)
    args = parser.parse_args()
    labels = cluster(args.indir, args.outpath)
    for label, id in labels.iteritems():
            print(label, id)
