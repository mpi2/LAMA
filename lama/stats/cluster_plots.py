
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
import matplotlib
from collections import OrderedDict
from os.path import isdir, split
from pathlib import Path
from typing import List, Union

# Force matplotlib to not use any Xwindows backend.
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import umap

from lama.common import LoadImage
from lama.stats.standard_stats.data_loaders import LineData

THRESH = 4.0
MAX_Z = 400  # If size extent is larger, then rescale to speed things up/reduce memory usage
PCA_N_COMPONENTS = 50

# t-sne parameters
TSNE_PARAMETERS = {
    'n_components': 2,
    'init': 'pca',
    'random_state': 7,
    'n_iter': 500,
    'perplexity': 4,
    'learning_rate': 100,
    'method': 'barnes_hut',
    'metric': 'correlation'
}


def _plot(data: pd.DataFrame, title: str, outpath):
    """
    Plt the results of the clustering
    """

    ax = sns.scatterplot(x='x', y='y', data=data)
    plt.title(title)

    i = 1

    id_map = []

    for spec_id, row in data.iterrows():
        ax.text(row['x'] + 0.08, row['y'], str(i), horizontalalignment='center', size='medium', color='black', weight='semibold')
        id_map.append([i, spec_id])
        i += 1

    plt.savefig(outpath)
    plt.close()

    df_ids = pd.DataFrame.from_records(id_map, columns=['id', 'specimen'], index='id')
    df_ids.to_csv(outpath.with_suffix('.csv'))




def umap_organs(data: pd.DataFrame, outpath: Path, title=''):


    embedding = umap.UMAP(n_neighbors=2,
                          min_dist=0.2,
                          metric='correlation').fit_transform(data)


    df = pd.DataFrame(embedding[:, 0:2], index=data.index, columns=['x', 'y'])

    _plot(df, title, outpath)


def pca_preprocess(data: pd.DataFrame, outpath,  title='') -> pd.DataFrame:

    # We try to reduce the feature numbers doen to a more managable level (~50)
    # However if we have less features than 50 (maybe during testing) set components to that level

    # If dataframe, extract the numpy array

    n_components = min(data.shape[1], data.shape[0])
    print(n_components)
    pca = PCA(n_components=n_components)

    try:
        pca_result = pca.fit_transform(data.values)
    except Exception:
        print(outpath.name, 'failed')
        return

    df = pd.DataFrame(pca_result[:, 0:2], index=data.index, columns=['x', 'y'])
    _plot(df, title, outpath)






