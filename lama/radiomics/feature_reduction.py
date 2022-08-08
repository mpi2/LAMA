import logging

from sklearn.ensemble import RandomForestClassifier
import numpy as np
from sklearn.feature_selection import SequentialFeatureSelector, SelectFromModel
from pathlib import Path
from itertools import product
import pandas as pd
def feat_select(X, k, score_method):
    m = SequentialFeatureSelector(RandomForestClassifier(), n_features_to_select=k, scoring=score_method)
    m.fit(X, X.index)
    return m.score(X, X.index)


def main():
    logging.info("Starting")
    X = pd.read_csv("Z:/jcsmr/ROLab/Experimental data/Radiomics/Workflow design and trial results/Kyle Drover analysis/220617_BQ_norm_stage_full/sub_normed_features.csv")

    X.set_index('Tumour_Model', inplace=True)

    X = X.select_dtypes(include=np.number)

    score_methods = ['roc_auc','precision', 'neg_log_loss', 'f1', 'recall', 'jaccard']
    k = [len(X.columns)-1, 1000, 500, 300, 200, 100, 50, 25, 10, 6, 4, 2, 1]
    k.reverse()
    all_combs = list(product(k, score_methods))
    print(all_combs)
    full_results = [None] * len(all_combs)

    for i, comb in enumerate(all_combs):
        logging.info(comb)
        print(comb)
        full_results[i] = [comb[0], comb[1], feat_select(X, comb[0], comb[1])]

    final_results = pd.DataFrame(full_results)
    final_results.to_csv("Z:/jcsmr/ROLab/Experimental data/Radiomics/Workflow design and trial results/Kyle Drover analysis/220617_BQ_norm_stage_full/feature_red_scores.csv")

    easy_m = SelectFromModel(estimator=RandomForestClassifier())
    easy_m.fit_transform(X, X.index)

    X.to_csv("Z:/jcsmr/ROLab/Experimental data/Radiomics/Workflow design and trial results/Kyle Drover analysis/220617_BQ_norm_stage_full/red_features.csv")

