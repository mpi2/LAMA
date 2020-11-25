import numpy as np
import pandas as pd
from lama import common
from typing import Union
from pathlib import Path


def folding_report(jac_array, label_map: np.ndarray, label_info: Union[pd.DataFrame, str, Path] = None, outdir=None):
    """
    Write out csv detailing the presence of folding per organ
    """
    if jac_array is None:
        return
    labels = np.unique(label_map)

    if not isinstance(label_info, pd.DataFrame):
        label_info = pd.read_csv(label_info, index_col=0)

    result = []
    for label in labels:
        # print(label)
        if label == 0:
            continue
        jac_label = jac_array[label_map == label]
        label_size = jac_label.size
        num_neg_vox = jac_label[jac_label < 0].size
        total_folding = np.sum(jac_label[jac_label < 0])
        result.append([label, label_size, num_neg_vox, total_folding])

    df = pd.DataFrame.from_records(result)
    df.columns = ['label', 'label_size', 'num_neg_voxels', 'summed_folding']
    df.set_index('label', drop=True, inplace=True)

    if label_info is not None:
        df = df.merge(label_info[['label_name']], left_index=True, right_index=True)

    if outdir:
        df.to_csv(outdir / common.FOLDING_FILE_NAME)
    else:
        return df