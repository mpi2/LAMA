from os.path import join, basename
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

from lama.common import LoadImage
from lama.stats.standard_stats.data_loaders import LineData

THRESH = 4.0
MAX_Z = 400  # If size extent is larger, then rescale to speed things up/reduce memory usage
PCA_N_COMPONENTS = 50

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


def tsne_on_raw_data(data: pd.DataFrame, out_dir: Path):
    """
    Given a list a mask-removed numpy arrays, cluster using t-sne

    Parameters
    ----------
    input_data
        Contains all the input data including masks and labels
    """

    # # Preprocess to reduce complexity and noise?
    # pca_data = _pca_preprocess(input_data.data)

    tsne = TSNE(**TSNE_PARAMETERS)

    trans_data = tsne.fit_transform(data.values)
    title = "t-SNE clustering on raw data"

    df = pd.DataFrame(trans_data, columns=['x', 'y'], index=list(data.index))

    sns.lmplot(x='x', y='y', data=df, fit_reg=False)

    # label_names_csv = join(split(outpath)[0], 'labels.csv')
    # with open(label_names_csv, 'w') as lh:
    #     for i in range(trans_data[0].size):
    #         plt.annotate(list(names.keys())[i], xy=(trans_data[0][i], trans_data[1][i]))
    #         lh.write('{},{}\n'.format(i, names[i]))

    # fig_text = '\n'.join([x for x in names.values()])
    # plt.gcf().text(1.2, 1, fig_text, fontsize=8)
    plt.title(title)
    out_path = out_dir / 'tsne_clustering.png'
    plt.savefig(out_path, bbox_inches='tight',dpi=100)
    plt.close()


def _pca_preprocess(data: Union[pd.DataFrame, np.ndarray]) -> pd.DataFrame:

    # We try to reduce the feature numbers doen to a more managable level (~50)
    # However if we have less features than 50 (maybe during testing) set components to that level

    # If dataframe, extract the numpy array
    try:
        data = data.values
    except AttributeError:
        pass

    n_components = min(data.shape[1], PCA_N_COMPONENTS)
    print(n_components)
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(data)
    print(f'Cumulative explained variation for 50 principal components: {np.sum(pca.explained_variance_ratio_)}')
    return pca_result



