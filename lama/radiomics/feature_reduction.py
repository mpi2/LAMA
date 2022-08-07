from sklearn.ensemble import RandomForestClassifier

from sklearn.feature_selection import SequentialFeatureSelector, SelectFromModel
from pathlib import Path
from itertools import product
import pandas as pd
def feat_select(X, k, score_method):
    m = SequentialFeatureSelector(RandomForestClassifier(), n_features_to_select=k, scoring=score_method)
    m.fit(X)
    return m.score(X)


def main():
    X = pd.read_csv("")
    score_methods = ['roc_auc','precision', 'neg_log_loss', 'f1', 'recall', 'jaccard']
    k = [len(X), 1000, 500, 300, 200, 100, 50, 25, 10, 6, 4, 2, 1]
    all_combs = product([k, score_methods])
    full_results = [None] * len(all_combs)

    for i, comb in enumerate(all_combs):
        full_results[i] = [comb[0], comb[1], feat_select(X, comb[0], comb[1])]

    final_results = pd.Dataframe(full_results)
    final_results.to_csv()

    easy_m = SelectFromModel(estimator=RandomForestClassifier())
    easy_m.fit_transform(X)

    X.to_csv()


