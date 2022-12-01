from logzero import logger as logging
import os
from catboost import CatBoostClassifier, Pool, sum_models, cv
import matplotlib.pyplot as plt
import time
import shap
#from mlxtend.feature_selection import SequentialFeatureSelector as SFS
import pickle
#from mlxtend.plotting import plot_sequential_feature_selection as plot_sfs

from sklearn.ensemble import RandomForestClassifier
import numpy as np
from sklearn.feature_selection import SequentialFeatureSelector, SelectFromModel
from sklearn.model_selection import LeaveOneOut
from pathlib import Path
from itertools import product
import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import accuracy_score, roc_auc_score, precision_score, f1_score, recall_score, matthews_corrcoef, make_scorer
from sklearn.pipeline import Pipeline
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import OneSidedSelection
from imblearn.pipeline import Pipeline
import statistics
import seaborn as sns
from collections import Counter

import torch
def correlation(dataset: pd.DataFrame, _dir: Path=None, threshold: float = 0.9, org=None):
    """
    identifies correlated features in  a

    Parameters
    ----------
    dataset: pandas dataframem
    """
    if org:
        _dir = _dir / str(org)
        os.makedirs(_dir, exist_ok=True)

    col_corr = set()  # Set of all the names of correlated columns
    corr_matrix = dataset.corr(method="spearman")

    logging.info("saving corr matrix at {}".format(_dir))
    fig, ax = plt.subplots(figsize=[50, 50])
    #cm = sns.diverging_palette(250, 15, s=100, as_cmap=True)
    sns.heatmap(corr_matrix, ax=ax,
                cbar_kws={'label': "Absolute value of Spearman's correlation"},
                square=True)
    ax.figure.axes[-1].yaxis.label.set_size(22)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    plt.tight_layout()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    plt.savefig(str(_dir) + "/corr_matix.png")
    plt.close()
    # do the removal
    for i in range(len(corr_matrix.columns)):
        for j in range(i):
            if abs(corr_matrix.iloc[i, j]) > threshold: # we are interested in absolute coeff value
                colname = corr_matrix.columns[i]  # getting the name of column
                col_corr.add(colname)
    return col_corr

def shap_feature_ranking(data, shap_values, columns=[]):
    """
    From stack-overflow, return columns
    """
    if not columns: columns = data.columns.tolist()  # If columns are not given, take all columns

    c_idxs = []
    for column in columns: c_idxs.append(
        data.columns.get_loc(column))  # Get column locations for desired columns in given dataframe
    if isinstance(shap_values, list):  # If shap values is a list of arrays (i.e., several classes)
        means = [np.abs(shap_values[class_][:, c_idxs]).mean(axis=0) for class_ in
                 range(len(shap_values))]  # Compute mean shap values per class
        shap_means = np.sum(np.column_stack(means), 1)  # Sum of shap values over all classes
    else:  # Else there is only one 2D array of shap values
        assert len(shap_values.shape) == 2, 'Expected two-dimensional shap values array.'
        shap_means = np.abs(shap_values).mean(axis=0)

    # Put into dataframe along with columns and sort by shap_means, reset index to get ranking
    df_ranking = pd.DataFrame({'feature': columns, 'mean_shap_value': shap_means}).sort_values(by='mean_shap_value',
                                                                                               ascending=False).reset_index(
        drop=True)
    df_ranking.index += 1
    return df_ranking



def shap_feat_select(X, m, _dir, cut_off: float=-1,n_feat_cutoff: float=None, org: int=None):
    """

    """
    #m = RandomForestClassifier(n_jobs=-1, n_estimators=100, verbose=0, oob_score=True)

    org_dir = _dir / str(org)

    os.makedirs(org_dir, exist_ok=True)

    cut_off_dir = org_dir / str(cut_off)
    os.makedirs(cut_off_dir, exist_ok=True)

    # print("plotting intrinsic RF rank")
    # importances = m.feature_importances_
    # indices = np.argsort(importances)
    # features = X.columns
    # plt.title('Feature Importances')
    # plt.figure(figsize=(15,200))
    # plt.rc('ytick', labelsize=6)
    # plt.barh(range(len(indices)), importances[indices], color='b', align='center')
    # plt.yticks(range(len(indices)), [features[i] for i in indices])
    # plt.xlabel('Relative Importance')
    # plt.tight_layout()
    # plt.savefig(str(cut_off_dir)+"/rf_rank.png")
    # plt.close()


    explainer = shap.KernelExplainer(m.predict, X, verbose=False)
    shap_values = explainer.shap_values(X)
    shap.summary_plot(shap_values, X, show=False, max_display=20)

    shap_importance = shap_feature_ranking(X, shap_values)

    # just a flag to have features in the first place

    if n_feat_cutoff:
        shap_importance = shap_importance[0:n_feat_cutoff]
        os.makedirs(_dir / str(n_feat_cutoff), exist_ok=True)
        X = X[shap_importance['feature']]

        plt.tight_layout()

        plt.savefig(str(_dir / str(n_feat_cutoff)) + "/shap_feat_rank_plot.png")

        plt.close()

    if cut_off >= 0:
        shap_importance =  shap_importance[shap_importance['mean_shap_value'] > cut_off]
        X = X[shap_importance['feature']]

        plt.tight_layout()

        plt.savefig(str(cut_off_dir) + "/shap_feat_rank_plot.png")

        plt.close()
    return X

def smote_oversampling(X, k: int=6, max_non_targets: int=150):
    # gets the ratio of target to baseline
    non_targets = Counter(X.index)['0']
    targets = Counter(X.index)['1']
    obs_ratio = targets / non_targets
    logging.info("Original ratio of targets : non-targets = {}".format(obs_ratio))
    if (non_targets > max_non_targets) & (obs_ratio < 0.2):
        required_ratio = 150 / non_targets
        logging.info("Undersampling to 150 targets to improve SHAP speed, followed by oversampling")
        steps = [('u', OneSidedSelection(sampling_strategy=required_ratio)),
                 ('o', SMOTE(n_jobs=-1, k_neighbors=k-1))]
        pipeline = Pipeline(steps=steps)
        x_train_std_os, y_train_os = pipeline.fit_resample(X, X.index)
    elif 0.9 <= obs_ratio <= 1.1:
        logging.info("dataset is relatively balanced, returning original data")
        x_train_std_os = X
    else:
        sm = SMOTE(n_jobs=-1, k_neighbors=k - 1)
        x_train_std_os, y_train_os = sm.fit_resample(X, X.index)

    x_train_std_os.set_index(y_train_os, inplace=True)
    logging.info("Rebalanced dataset: {}".format(Counter(x_train_std_os.index)))

    return x_train_std_os


def main(X, org, rad_file_path, batch_test = None):

    logging.info("Doing org: {}".format(org))

    print("Doing org: {}".format(org))
    logging.info("Starting")

    #X = X[X['org']== org]


    if org:
        X['HPE']  = X['HPE'].map({'normal': 0, 'abnormal': 1}).astype(int)
        X.set_index('HPE', inplace=True)
        #X = X[X['background'] == 'C3HHEH']
        #X['genotype'] = X['genotype'].map({'WT': 0, 'HET': 1}).astype(int)
        #X.set_index('genotype', inplace=True)

    elif batch_test:
        X = X[(X['Age'] == 'D14') & (X['Tumour_Model'] == '4T1R')]
        X['Exp'] = X['Exp'].map({'MPTLVo4': 0, 'MPTLVo7': 1})
        X.set_index('Exp', inplace=True)
        X.drop(['Date', 'Animal_No.'], axis=1, inplace=True)

    else:
        X['Tumour_Model'] = X['Tumour_Model'].map({'4T1R': 0, 'CT26R': 1}).astype(int)
        X.set_index('Tumour_Model', inplace=True)
        X.drop(['Date', 'Animal_No.'], axis=1, inplace=True)



    X = X.select_dtypes(include=np.number)

    # lets remove correlated variables

    corr_feats = correlation(X, rad_file_path.parent, 0.9, org=org)

    logging.info('{}: {}'.format("Number of features removed due to correlation", len(set(corr_feats))))


    X.drop(corr_feats, axis=1, inplace=True)


    # clone X for final test
    X_to_test = X

    #shap_cut_offs = [0.005, 0.01, 0.02]

    org_dir = _dir=rad_file_path.parent / str(org)


    #parameters = {'num_trees': list(range(50, 1000, 100))}



    #X=shap_feat_select(X)

    # balancing clsses via SMOTE
    logging.info("oversampling via smote")
    n_test = X[X.index == 1].shape[0]

    X = smote_oversampling(X, n_test) if n_test < 5 else smote_oversampling(X)

    logging.info("fitting model to training data")
    m = CatBoostClassifier(iterations=1000, task_type='GPU', verbose=250, train_dir=org_dir)
    m.fit(X, X.index.to_numpy())
    logging.info("doing feature selection using SHAP")

    if org:
        shap_cut_offs = list(np.arange(0.000, 0.025, 0.005))
        full_X = [shap_feat_select(X, m, _dir=rad_file_path.parent, cut_off=cut_off, org=org) for cut_off in shap_cut_offs]
    else:
        n_feats = list(np.arange(20, 21, 1))
        full_X = [shap_feat_select(X, m, _dir=rad_file_path.parent, n_feat_cutoff=n, org=org) for n in n_feats]


    n_feats = [X.shape[1] for X in full_X]

    logging.info("n_feats: {}".format(n_feats))


    for i, x in enumerate(full_X):
        # make different models for different feature nums
        models = []
        model_dir = org_dir / str(x.shape[1])
        os.makedirs(model_dir, exist_ok=True)

        # sample 20 different train-test partitions (train size of 0.2) and create an average model
        for j in range(20):
            train_dir = model_dir / str(j)
            os.makedirs(train_dir, exist_ok=True)
            all_x = Pool(data=x, label=x.index.to_numpy())

            m = CatBoostClassifier(iterations=1000, task_type="CPU", loss_function='Logloss', train_dir=train_dir,
                                   custom_loss=['AUC', 'Accuracy', 'Precision', 'F1', 'Recall'],
                                   verbose=250)

            # optimise via grid search

            params = {
                'depth': [4, 6, 10],
                'l2_leaf_reg': [3, 5, 7],
            }

            m.grid_search(params, all_x, cv=20)

            logging.info("grid search: Number of trees {}, best_scores {}".format(m.tree_count_, m.get_best_score()))

            # now train with optimised parameters on split
            x_train, x_test, y_train, y_test = train_test_split(x, x.index.to_numpy(), test_size=0.20)

            train_pool = Pool(data=x_train, label=y_train)
            validation_pool = Pool(data=x_test, label=y_test)

            m.fit(train_pool, eval_set=validation_pool, verbose=False)

            logging.info("Eval CPU: Number of trees {}, best_scores {}".format(m.tree_count_, m.get_best_score()))

            # perform 30 fold cross-validation

            cv_data = cv(params=m.get_params(),
                         pool=train_pool,
                         fold_count=30,
                         shuffle=True,
                         partition_random_seed=0,
                         stratified=False,
                         verbose=True,
                         plot=False,
                         as_pandas=True,
                         return_models=False)

            cv_filename = str(rad_file_path.parent) + "/" + str(len(set(corr_feats))) + "_" + str(j) + ".csv"

            logging.info("saving cv results to {}".format(cv_filename))
            cv_data.to_csv(cv_filename)

            # tests GPU training
            #TODO: remove if useless
            m2 = CatBoostClassifier(iterations=1000, task_type="GPU",  train_dir=str(train_dir),
                                    custom_loss=['Accuracy', 'Precision', 'F1', 'Recall'],
                                    verbose=250)
            m2.fit(train_pool, eval_set=validation_pool, verbose=False)
            logging.info("Eval GPU: Number of trees {}, best_scores {}".format(m2.tree_count_, m2.get_best_score()))

            logging.info("Saving models")
            m_filename = str(rad_file_path.parent) + "/" + str(org) + "/CPU_" + str(x.shape[1]) + "_" + str(j) + ".cbm"
            m2_filename = str(rad_file_path.parent) + "/" + str(org) + "/GPU_" + + str(x.shape[1]) + "_" + str(j) + ".cbm"

            m.save_model(m_filename)
            m2.save_model(m2_filename)

            models.append(m)
            models.append(m2)

        logging.info("Combining model predictions into one mega model")

        m_avg = sum_models(models, weights=[1.0 / len(models)] * len(models))

        avrg_filename = str(rad_file_path.parent) + "/" + str(org) + "/mega_model" + str(x.shape[1]) + ".cbm"

        m_avg.save_model(avrg_filename)

        logging.info("Mega_Model: Number of trees {}, best_scores {}".format(m_avg.tree_count_, m_avg.get_best_score()))





