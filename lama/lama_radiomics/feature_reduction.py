from logzero import logger as logging
import os

import matplotlib.pyplot as plt
import time
import shap
from mlxtend.feature_selection import SequentialFeatureSelector as SFS
import pickle
from mlxtend.plotting import plot_sequential_feature_selection as plot_sfs

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
import statistics

def correlation(dataset: pd.DataFrame, threshold: float = 0.9):
    """
    identifies correlated features in  a

    Parameters
    ----------
    dataset: pandas dataframem
    """
    col_corr = set()  # Set of all the names of correlated columns
    corr_matrix = dataset.corr(method="spearman")
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



def shap_feat_select(X, _dir, cut_off: float=-1, org: int=None):
    """

    """
    m = RandomForestClassifier(n_jobs=-1, n_estimators=100, verbose=0, oob_score=True)
    print("fitting model to training data")
    m.fit(X, X.index)
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
    if cut_off >= 0:
        shap_importance =  shap_importance[shap_importance['mean_shap_value'] > cut_off]

    X = X[shap_importance['feature']]

    plt.tight_layout()

    plt.savefig(str(cut_off_dir) + "/shap_feat_rank_plot.png")

    plt.close()
    return X



def smote_oversampling(X, k: int=6):
    sm = SMOTE(n_jobs=-1, k_neighbors=k-1)
    x_train_std_os, y_train_os = sm.fit_resample(X, X.index)
    x_train_std_os.set_index(y_train_os, inplace=True)
    return  x_train_std_os




def main(X, org, rad_file_path):

    logging.info("Doing org: {}".format(org))

    print("Doing org: {}".format(org))
    logging.info("Starting")

    #X = X[X['org']== org]

    if org:
        X['HPE']  = X['HPE'].map({'normal': 0, 'abnormal': 1}).astype(int)
        X.set_index('HPE', inplace=True)
    else:
        X['Tumour_Model'] = X['Tumour_Model'].map({'4T1R': 0, 'CT26R': 1}).astype(int)
        X.set_index('Tumour_Model', inplace=True)
        X.drop(['Date', 'Animal_No.'], axis=1, inplace=True)


    #X = pd.read_csv("Z:/jcsmr/ROLab/Experimental data/Radiomics/Workflow design and trial results/Kyle Drover analysis/220617_BQ_norm_stage_full/sub_normed_features.csv")

    #X['Tumour_Model'] = X['Tumour_Model'].map({'4T1R': 0, 'CT26R': 1}).astype(int)
    #X.set_index('Tumour_Model', inplace=True)
    #X.drop(['Date', 'scanID', 'Animal_No.'], axis=1, inplace=True)

    X = X.select_dtypes(include=np.number)


    # lets remove correlated variables

    corr_feats = correlation(X, 0.9)

    logging.info('{}: {}'.format("Number of features removed due to correlation", len(set(corr_feats))))


    X.drop(corr_feats, axis=1, inplace=True)



    #shap_cut_offs = [0.005, 0.01, 0.02]

    org_dir = _dir=rad_file_path.parent / str(org)

    parameters = {'n_estimators': list(range(50, 1000, 100))}



    #X=shap_feat_select(X)

    # balancing classes via SMOTE
    logging.info("oversampling via smote")
    n_test = X[X.index == 1].shape[0]

    X = smote_oversampling(X, n_test) if n_test < 5 else smote_oversampling(X)

    logging.info("doing feature selection using SHAP")

    shap_cut_offs = list(np.arange(0.005, 0.02, 0.005))

    full_X = [shap_feat_select(X, _dir=rad_file_path.parent, cut_off=cut_off, org=org) for cut_off in shap_cut_offs]


    n_feats = [X.shape[1] for X in full_X]

    logging.info("n_feats: {}".format(n_feats))

    scoring = {'AUC': 'roc_auc',
               'Accuracy': make_scorer(accuracy_score),
               'Precision': make_scorer(precision_score),
               'F1': make_scorer(f1_score),
               'Recall': make_scorer(recall_score),
               'MCC': make_scorer(matthews_corrcoef)}



    results = [None] * len(shap_cut_offs)

    results_v2 = [None] * len(shap_cut_offs)

    X_axises = [None] * len(shap_cut_offs)

    for i, x in enumerate(full_X):

        gs = GridSearchCV(RandomForestClassifier(oob_score=True),
                          param_grid=parameters,
                          n_jobs=-1,
                          scoring=scoring,
                          cv=20,
                          refit='AUC', verbose=1)

        gs.fit(x, x.index)
        results[i] = gs.cv_results_


        X_axises[i] = np.array(results[i]['param_n_estimators'].data, dtype=float)

    parameters = {'n_estimators': [175]}

    for i, x in enumerate(full_X):

        gs = GridSearchCV(RandomForestClassifier(oob_score=True),
                          param_grid=parameters,
                          n_jobs=-1,
                          scoring=scoring,
                          cv=20,
                          refit='AUC', verbose=1)

        gs.fit(x, x.index)
        results_v2[i] = gs.cv_results_



    colours =  ['k', 'b', 'c', 'm', 'r', 'y', 'g']

    n_trees_lst = []
    # so I want to make a different plot per metric
    for scorer in scoring.keys():

        plt.figure(figsize=(20, 20))
        plt.title("GridSearchCV evaluating using multiple scorers simultaneously",
                  fontsize=16)

        plt.xlabel("Number of Trees")
        plt.ylabel("Score")
        plt.grid()

        ax = plt.axes()
        ax.set_xlim(0, 1000)
        ax.set_ylim(0, 1.05)

        for i, (x_axis, colour) in enumerate(zip(X_axises, colours)):

            sample_score_mean = results[i]['mean_test_%s' % (scorer)]
            sample_score_std = results[i]['std_test_%s' % (scorer)]

            ax.fill_between(x_axis, sample_score_mean - sample_score_std,
                            sample_score_mean + sample_score_std,
                            alpha=0.1, color=colour)

            ax.plot(x_axis, sample_score_mean, '--', color=colour,
                    alpha=1,
                    label="%s using a SHAP cut-off of %s" % (scorer, shap_cut_offs[i]))
            best_index = np.nonzero(results[i]['rank_test_%s' % scorer] == 1)[0][0]
            best_score = results[i]['mean_test_%s' % scorer][best_index]

            # Plot a dotted vertical line at the best score for that scorer marked by x
            ax.plot([x_axis[best_index], ] * 2, [0, best_score],
                    linestyle='-.', color=colour, marker='x', markeredgewidth=3, ms=8)

            ax.annotate("%0.2f" % best_score,
                        (x_axis[best_index], best_score + 0.005))
            n_trees_lst.append(x_axis[best_index])


        plt.legend(loc="best")
        if org:
            plt.savefig(str(org_dir) + "/" + str(scorer) + "_curve.png")
        else:
            plt.savefig(str(rad_file_path.parent) + "/" + str(scorer) + "_curve.png")
        plt.close()

    best_ntrees = statistics.mode(n_trees_lst)

    print("best_ntrees: ", best_ntrees)
    x_axis = n_feats

    results = pd.DataFrame(results_v2).drop(columns='params').applymap(lambda x: float(x))

    plt.figure(figsize=(20, 20))
    plt.title("GridSearchCV evaluating using multiple scorers simultaneously",
              fontsize=16)

    plt.xlabel("Number of Features")
    plt.ylabel("Score")
    plt.grid()
    ax = plt.axes()
    ax.set_xlim(0, 150)
    ax.set_ylim(0, 1.05)


    best_cutoff_lst = []
    for scorer, colour in zip(scoring, colours):

        sample_score_mean = results['mean_test_%s' % (scorer)]

        sample_score_std = results['std_test_%s' % (scorer)]

        ax.fill_between(x_axis, sample_score_mean - sample_score_std,
                            sample_score_mean + sample_score_std,
                            alpha=0.1, color=colour)

        ax.plot(x_axis, sample_score_mean, '--', color=colour,
                    alpha=1,
                    label="%s" % (scorer))
        best_index =  np.argmax(sample_score_mean)

        best_score = results['mean_test_%s' % (scorer)][best_index]

        logging.info("best score: {} {}".format(scorer, best_score))


        # Plot a dotted vertical line at the best score for that scorer marked by x
        ax.plot([x_axis[best_index], ] * 2, [0, best_score],
                linestyle='-.', color=colour, marker='x', markeredgewidth=3, ms=8)

        ax.annotate("%0.2f" % best_score,
                    (x_axis[best_index], best_score + 0.005))


        best_cutoff_lst.append(x_axis[best_index])

    plt.legend(loc="best")
    if org:
        plt.savefig(str(org_dir) + "/feat_test_curve.png")
    else:
        plt.savefig(str(rad_file_path.parent) + "/feat_test_curve.png")
    plt.close()

    best_cut_off = statistics.mode(best_cutoff_lst)
    print(best_cut_off)

    print(n_feats.index(best_cut_off))

    best_X = full_X[n_feats.index(best_cut_off)]

    logging.info("Splitting data into 80-20 split")
    x_train, x_test, y_train, y_test = train_test_split(best_X, best_X.index, test_size=0.2, random_state=0)

    model = RandomForestClassifier(n_estimators=int(best_ntrees))

    model.fit(x_train, y_train)
    model_file_name = str(org_dir / "finalised_model.sav")

    pickle.dump(model, open(model_file_name, 'wb'))

    logging.info("final model score: {}". format(model.score(x_test, y_test)))

