#!/usr/bin/env python

import SimpleITK as sitk
from os.path import join, basename
from sklearn.manifold import TSNE
import numpy as np
import pandas as pd
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
from os.path import isdir, split
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common
from common import LoadImage
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



def cluster_from_array(imgs, ids, outpath, groups=None, label_map=None):
    """
    Given a list a mask-removed numpy arrays, cluster using t-sne

    Parameters
    ----------
    imgs: list
        list of arrays of masked image data
    ids: list
        ids associated with each img array
    groups: pandas.DataFrame
        column for id_ and group
        id_ is the basename of the image (no extension) group is for the hue in the plot (eg WT, baseline)
    masked_labels:
        non-mask regions of labels map
    """
    names = OrderedDict()
    for i, n in enumerate(ids):
        names[i] = n
    return _make_plot(imgs, names, outpath, groups, label_map)


def cluster_form_directory(indir, outpath):
    if isdir(indir):
        names = OrderedDict()
        paths = common.get_file_paths(indir)
    else:
        paths = common.get_inputs_from_file_list(indir)

    imgs = []
    first_img = LoadImage(paths[0]).array
    if first_img.shape[0] > MAX_Z:
        # determine scaling factor
        scale_factor = int(first_img.shape[0] / MAX_Z)
        if scale_factor < 2:
            scale_factor = 2
    else:
        scale_factor = False

    for i, vp in enumerate(paths):
        print(vp)
        bn = basename(vp)
        names[i] = bn
        img = LoadImage(vp).itkimg
        if scale_factor:
            img = sitk.BinShrink(img, (scale_factor, scale_factor, scale_factor))
        arr = sitk.GetArrayFromImage(img)
        try:
            arr[np.where((arr > -THRESH) & (arr < THRESH))] = 0
        except ValueError:  # No values less than 4
            continue
        imgs.append(arr.ravel())

    return _make_plot(imgs, names, outpath)  # Return the image names so they can be added to the log (should just put them in a legend on the figure instead)


def _make_plot(imgs, names, outpath, groups=None, label_map=None):
    """

    Parameters
    ----------
    imgs: list
        List a ndarrays of the data to cluster on. One speciemn per element
    names:

    outpath: str
        Path for the plot output
    groups: pandas.DataFrame
        group membership for each specimen

    Returns
    -------

    """

    import seaborn as sns
    import pandas as pd

    tsne = TSNE(**TSNE_PARAMETERS)

    if label_map is not None:
        thresh = 1
        label_hit_counts = []
        for img in imgs:
            img_result = []
            for label in range(1, label_map.max() + 1):
                # TODO: I think this could be sped up
                mean = np.mean(np.abs(img[label_map==label]))
                if np.isnan(mean):
                    mean = 0
                img_result.append(mean)
            label_hit_counts.append(np.array(img_result))
        trans_data = tsne.fit_transform(label_hit_counts).T
        title = "t-SNE clustering on mean label hit score"

    else:
        trans_data = tsne.fit_transform(imgs).T
        title = "t-SNE clustering on z-score"

    try:
        df = pd.DataFrame(trans_data.T, columns=['x', 'y'], index=names.values())
        if groups is not None:
            df['groups'] = groups['group'].values

    except Exception as e:
        pass # testing

    if groups is not None:
        try:
            sns.lmplot(x='x', y='y', data=df, fit_reg=False, hue='groups')
        except Exception as e:
            print(e)

    else:
        sns.lmplot(x='x', y='y', data=df, fit_reg=False)

    label_names_csv = join(split(outpath)[0], 'labels.csv')
    with open(label_names_csv, 'w') as lh:
        for i in range(trans_data[0].size):
            plt.annotate(names.keys()[i], xy=(trans_data[0][i], trans_data[1][i]))
            lh.write('{},{}\n'.format(i, names[i]))

    # fig_text = '\n'.join([x for x in names.values()])
    # plt.gcf().text(1.2, 1, fig_text, fontsize=8)
    plt.title(title)

    plt.savefig(outpath, bbox_inches='tight',dpi=100)
    plt.close()

    return names


def load_from_csv(in_path, outdir):
    import SimpleITK as sitk
    from scipy.stats import zmap
    import re

    df = pd.read_csv(in_path)
    arrays = []
    ids = []
    groups = []
    for index, row in df.iterrows():
        path = row['path']

        m = re.search('volumes/(.+?)_download', path)
        # m = re.search('jacobians/(.+?)_rec', path, re.IGNORECASE)
        ids.append(m.group(1))
        genotype = row['genotype']
        groups.append(genotype)
        arr = LoadImage(path).array.ravel()
        arrays.append(arr)
    # make dataframe for cluster function
    df_groups = pd.DataFrame.from_dict(dict(id_=ids, group=groups))
    # Now zmap
    arrays = np.vstack(arrays)
    zscores = []
    for a in arrays:
        z = zmap(a, arrays)
        zscores.append(z)
    zscores = np.vstack(zscores)
    outpath = join(outdir, 'cluster.png')
    cluster_from_array(zscores, ids, outpath, df_groups)

#
if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("Create t-sne clustering plot of images")
    parser.add_argument('-i', dest='in_csv', help='csv file with paths and group membership', required=True)
    parser.add_argument('-o', dest='outdir', help='outdir to stick the plot(s)', required=True)

    args = parser.parse_args()
    labels = load_from_csv(args.in_csv, args.outdir)

