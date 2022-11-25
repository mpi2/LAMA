from logzero import logger as logging
import os
from catboost import CatBoostClassifier, Pool
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
import statistics
import seaborn as sns
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



def smote_oversampling(X, k: int=6):
    sm = SMOTE(n_jobs=-1, k_neighbors=k-1)
    x_train_std_os, y_train_os = sm.fit_resample(X, X.index)
    x_train_std_os.set_index(y_train_os, inplace=True)
    return  x_train_std_os




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
    m = CatBoostClassifier(iterations=1000, task_type='GPU', verbose=100)
    m.fit(X, X.index.to_numpy())
    logging.info("doing feature selection using SHAP")
    #
    if org:
        shap_cut_offs = list(np.arange(0.000, 0.025, 0.005))
        full_X = [shap_feat_select(X, m, _dir=rad_file_path.parent, cut_off=cut_off, org=org) for cut_off in shap_cut_offs]
    else:
        n_feats = list(np.arange(20, 21, 1))
        full_X = [shap_feat_select(X, m, _dir=rad_file_path.parent, n_feat_cutoff=n, org=org) for n in n_feats]


    n_feats = [X.shape[1] for X in full_X]

    logging.info("n_feats: {}".format(n_feats))

    scoring = {'AUC': 'roc_auc',
               'Accuracy': make_scorer(accuracy_score),
               'Precision': make_scorer(precision_score),
               'F1': make_scorer(f1_score),
               'Recall': make_scorer(recall_score),
               'MCC': make_scorer(matthews_corrcoef)}



    results = [None] * len(n_feats)

    results_v2 = [None] * len(n_feats)

    X_axises = [None] * len(n_feats)

    parameters = {'iterations': list(range(200, 1000, 400))}
    for i, x in enumerate(full_X):
        for i in range(20):
            all_x = Pool(data=x, label=x.index.to_numpy())

            m = CatBoostClassifier(iterations=1000, task_type="GPU", custom_loss=['AUC', 'Accuracy','Precision', 'F1', 'Recall'],
                                   verbose=100)

            params = {
                'iterations': [600, 1000, 1400],
                'depth': [4, 6, 10],
                'l2_leaf_reg': [1, 3, 5, 7],
                }

            m.grid_search(params, all_x, cv=5)

            print(m.get_best_score())

            x_train, x_test, y_train, y_test = train_test_split(x, x.index.to_numpy(), test_size=0.20)

            train_pool = Pool(data=x_train, label=y_train)
            validation_pool = Pool(data=x_test, label=y_test)

            m.fit(train_pool, eval_set=validation_pool, verbose=False)
            print(m.get_best_score())
            print(m.get_evals_result())

        return True

        from catboost import cv




        cv_data = cv(params=params,
                     pool=train_pool,
                     fold_count=2,
                     shuffle=True,
                     partition_random_seed=0,
                     plot=False,
                     stratified=False,
                     verbose=True)
        print('cv_data', cv_data)


        gs = GridSearchCV(CatBoostClassifier(task_type="GPU",  verbose=False),
                                            param_grid=parameters,
                                            n_jobs=2,
                                            scoring=scoring,
                                            cv=LeaveOneOut(),
                                            refit='Accuracy', verbose=0)

        gs.fit(x, x.index.to_numpy())
        print(gs.best_params_)
        print(gs.best_score_)

        # split data???

        #gs = GridSearchCV(CatBoostClassifier(task_type="CPU"),
        #                  param_grid=parameters,
        #                  n_jobs=-1,
        #                  scoring=scoring,
        #                  cv=10,
        #                  refit='Accuracy', verbose=0)
        #gs =  CatBoostClassifier(iterations=1000, task_type="CPU")
        #gs_result = gs.grid_search(parameters, x, x.index.to_numpy(), plot=True)


        #gs.fit(x, x.index.to_numpy())
        #results[i] = gs.cv_results_
        #print("iterations results: ", gs.cv_results_)

        X_axises[i] = np.array(results[i]['param_iterations'].data, dtype=float)

    parameters = {'iterations': [1000]}

    for i, x in enumerate(full_X):

        gs = GridSearchCV(CatBoostClassifier(task_type="GPU"),
                          param_grid=parameters,
                          n_jobs=-1,
                          scoring=scoring,
                          cv=5,
                          refit='AUC', verbose=0)

        gs.fit(x, x.index.to_numpy())
        print(gs.cv_results_)
        print(gs.best_score_)
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
            print("iterations sample score", sample_score_mean)
            sample_score_std = results[i]['std_test_%s' % (scorer)]

            ax.fill_between(x_axis, sample_score_mean - sample_score_std,
                            sample_score_mean + sample_score_std,
                            alpha=0.1, color=colour)

            ax.plot(x_axis, sample_score_mean, '--', color=colour,
                    alpha=1,
                    label="%s number of feats %s" % (scorer, n_feats[i]))
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
    print("x_axis", x_axis)
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

        print("sample_score_mean", sample_score_mean)
        print(best_index)

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

    logging.info("best num of feats: {}".format(n_feats[n_feats.index(best_cut_off)]))

    best_X = full_X[n_feats.index(best_cut_off)]

    logging.info("doing final SHAP")
    X_to_test= shap_feat_select(X_to_test, _dir=rad_file_path.parent, n_feat_cutoff=n_feats[n_feats.index(best_cut_off)],
                               org=org)

    logging.info("Splitting data into 80-20 split")
    x_train, x_test, y_train, y_test = train_test_split(X_to_test, X_to_test.index, test_size=0.2, random_state=0)


    fig, ax = plt.subplots(figsize=[50, 50])
    sns.pairplot(x_train)

    ax.figure.axes[-1].yaxis.label.set_size(22)

    plt.tight_layout()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    plt.savefig(str(rad_file_path.parent) + "/pair_plot.png")
    plt.close()

    for x in range(20):
        model = CatBoostClassifier(iterations=best_ntrees, random_state=x)
        model.fit(x_train, y_train)
        model_file_name = str(org_dir / "finalised_model.sav")
        pickle.dump(model, open(model_file_name, 'wb'))
        #logging.info("final model score: {}". format(model.score(x_test, y_test)))
        logging.info("final ROC AUC score: {}".format(accuracy_score(y_test, model.predict(x_test))))
        #logging.info("decision path: {}".format(model.decision_path(x_test)))
        #logging.info("leaf indices: {}".format(model.apply(x_test)))

