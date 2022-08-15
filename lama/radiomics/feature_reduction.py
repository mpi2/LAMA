import logging
import matplotlib.pyplot as plt

from mlxtend.feature_selection import SequentialFeatureSelector as SFS

from mlxtend.plotting import plot_sequential_feature_selection as plot_sfs

from sklearn.ensemble import RandomForestClassifier
import numpy as np
from sklearn.feature_selection import SequentialFeatureSelector, SelectFromModel
from pathlib import Path
from itertools import product
import pandas as pd
from sklearn.model_selection import train_test_split
def feat_select(X, score_method):
    print("doing ", score_method)
    m = RandomForestClassifier()
    print("fitting model")
    m.fit(X, X.index)
    print("doing sfs")
    sfs1 = SFS(m, k_features=(1, 40), floating=True, forward=True, scoring=score_method, n_jobs=-1, verbose=2)

    sfs1.fit(X,x.index)

    print("getting metric dict")
    metric_dict = sfs1.get_metric_dict(confidence_interval=0.95)
    print(sfs1.k_feature_names_)

    print("plotting")
    plot_sfs(metric_dict)

    plt.savefig("Z:/jcsmr/ROLab/Experimental data/Radiomics/Workflow design and trial results/Kyle Drover analysis/220617_BQ_norm_stage_full/"+str(score_method)+".png")

    return sfs1.subsets_


def main():
    logging.info("Starting")
    X = pd.read_csv("Z:/jcsmr/ROLab/Experimental data/Radiomics/Workflow design and trial results/Kyle Drover analysis/220617_BQ_norm_stage_full/sub_normed_features.csv")

    X['Tumour_Model'] = X['Tumour_Model'].map({'4T1R': 0, 'CT26R': 1}).astype(int)
    X.set_index('Tumour_Model', inplace=True)
    X.drop(['Date'], axis=1, inplace=True)
    X.drop(['scanID'], axis=1, inplace=True)
    X.drop(['Animal_No.'], axis=1, inplace=True)
    print(X)
    X = X.select_dtypes(include=np.number)


    score_methods = ['f1', 'accuracy', 'roc_auc','precision', 'recall']
    k = [len(X.columns)-1, 1000, 500, 300, 200, 100, 50, 25, 10, 6, 4, 2, 1]
    k.reverse()
    all_combs = list(product(k, score_methods))
    full_results = [None] * len(all_combs)

    for i, meth in enumerate(score_methods):
        full_results[i] = [meth, feat_select(X, meth)]
    #for i, comb in enumerate(all_combs):
    #    logging.info(comb)
    #    print(comb)
    #    full_results[i] = [comb[0], comb[1], feat_select(X, comb[0], comb[1])]

    final_results = pd.DataFrame(full_results)
    final_results.to_csv("Z:/jcsmr/ROLab/Experimental data/Radiomics/Workflow design and trial results/Kyle Drover analysis/220617_BQ_norm_stage_full/feature_red_scores.csv")

    easy_m = SelectFromModel(estimator=RandomForestClassifier())
    easy_results = easy_m.fit_transform(X, X.index)

    easy_results.to_csv("Z:/jcsmr/ROLab/Experimental data/Radiomics/Workflow design and trial results/Kyle Drover analysis/220617_BQ_norm_stage_full/red_features.csv")


