
# print(m.get_evals_result())




# scoring = {'AUC': 'roc_auc',
    #           'Accuracy': make_scorer(accuracy_score),
    #           'Precision': make_scorer(precision_score),
    #           'F1': make_scorer(f1_score),
    #           'Recall': make_scorer(recall_score),
    #           'MCC': make_scorer(matthews_corrcoef)}

    # results = [None] * len(n_feats)
    # results_v2 = [None] * len(n_feats)
    # X_axises = [None] * len(n_feats)
    # parameters = {'iterations': list(range(200, 1000, 400))}

# split data???

# gs = GridSearchCV(CatBoostClassifier(task_type="CPU"),
#                  param_grid=parameters,
#                  n_jobs=-1,
#                  scoring=scoring,
#                  cv=10,
#                  refit='Accuracy', verbose=0)
# gs =  CatBoostClassifier(iterations=1000, task_type="CPU")
# gs_result = gs.grid_search(parameters, x, x.index.to_numpy(), plot=True)


# gs.fit(x, x.index.to_numpy())
# results[i] = gs.cv_results_
# print("iterations results: ", gs.cv_results_)

# X_axises[i] = np.array(results[i]['param_iterations'].data, dtype=float)
#
# parameters = {'iterations': [1000]}
#
# for i, x in enumerate(full_X):
#
#     gs = GridSearchCV(CatBoostClassifier(task_type="GPU"),
#                       param_grid=parameters,
#                       n_jobs=-1,
#                       scoring=scoring,
#                       cv=5,
#                       refit='AUC', verbose=0)
#
#     gs.fit(x, x.index.to_numpy())
#     print(gs.cv_results_)
#     print(gs.best_score_)
#     results_v2[i] = gs.cv_results_
#
#
#
# colours =  ['k', 'b', 'c', 'm', 'r', 'y', 'g']
#
# n_trees_lst = []
# # so I want to make a different plot per metric
# for scorer in scoring.keys():
#
#     plt.figure(figsize=(20, 20))
#     plt.title("GridSearchCV evaluating using multiple scorers simultaneously",
#               fontsize=16)
#
#     plt.xlabel("Number of Trees")
#     plt.ylabel("Score")
#     plt.grid()
#
#     ax = plt.axes()
#     ax.set_xlim(0, 1000)
#     ax.set_ylim(0, 1.05)
#
#     for i, (x_axis, colour) in enumerate(zip(X_axises, colours)):
#
#         sample_score_mean = results[i]['mean_test_%s' % (scorer)]
#         print("iterations sample score", sample_score_mean)
#         sample_score_std = results[i]['std_test_%s' % (scorer)]
#
#         ax.fill_between(x_axis, sample_score_mean - sample_score_std,
#                         sample_score_mean + sample_score_std,
#                         alpha=0.1, color=colour)
#
#         ax.plot(x_axis, sample_score_mean, '--', color=colour,
#                 alpha=1,
#                 label="%s number of feats %s" % (scorer, n_feats[i]))
#         best_index = np.nonzero(results[i]['rank_test_%s' % scorer] == 1)[0][0]
#         best_score = results[i]['mean_test_%s' % scorer][best_index]
#
#         # Plot a dotted vertical line at the best score for that scorer marked by x
#         ax.plot([x_axis[best_index], ] * 2, [0, best_score],
#                 linestyle='-.', color=colour, marker='x', markeredgewidth=3, ms=8)
#
#         ax.annotate("%0.2f" % best_score,
#                     (x_axis[best_index], best_score + 0.005))
#         n_trees_lst.append(x_axis[best_index])
#
#
#     plt.legend(loc="best")
#     if org:
#         plt.savefig(str(org_dir) + "/" + str(scorer) + "_curve.png")
#     else:
#         plt.savefig(str(rad_file_path.parent) + "/" + str(scorer) + "_curve.png")
#     plt.close()
#
# best_ntrees = statistics.mode(n_trees_lst)
#
# print("best_ntrees: ", best_ntrees)
# x_axis = n_feats
#
# results = pd.DataFrame(results_v2).drop(columns='params').applymap(lambda x: float(x))
#
# plt.figure(figsize=(20, 20))
# plt.title("GridSearchCV evaluating using multiple scorers simultaneously",
#           fontsize=16)
#
# plt.xlabel("Number of Features")
# plt.ylabel("Score")
# plt.grid()
# ax = plt.axes()
# ax.set_xlim(0, 150)
# ax.set_ylim(0, 1.05)
#
#
# best_cutoff_lst = []
# print("x_axis", x_axis)
# for scorer, colour in zip(scoring, colours):
#
#     sample_score_mean = results['mean_test_%s' % (scorer)]
#
#     sample_score_std = results['std_test_%s' % (scorer)]
#
#     ax.fill_between(x_axis, sample_score_mean - sample_score_std,
#                     sample_score_mean + sample_score_std,
#                     alpha=0.1, color=colour)
#
#     ax.plot(x_axis, sample_score_mean, '--', color=colour,
#             alpha=1,
#             label="%s" % (scorer))
#     best_index =  np.argmax(sample_score_mean)
#
#     print("sample_score_mean", sample_score_mean)
#     print(best_index)
#
#     best_score = results['mean_test_%s' % (scorer)][best_index]
#
#     logging.info("best score: {} {}".format(scorer, best_score))
#
#
#     # Plot a dotted vertical line at the best score for that scorer marked by x
#     ax.plot([x_axis[best_index], ] * 2, [0, best_score],
#             linestyle='-.', color=colour, marker='x', markeredgewidth=3, ms=8)
#
#     ax.annotate("%0.2f" % best_score,
#                 (x_axis[best_index], best_score + 0.005))
#
#
#     best_cutoff_lst.append(x_axis[best_index])
#
# plt.legend(loc="best")
# if org:
#     plt.savefig(str(org_dir) + "/feat_test_curve.png")
# else:
#     plt.savefig(str(rad_file_path.parent) + "/feat_test_curve.png")
# plt.close()
#
# best_cut_off = statistics.mode(best_cutoff_lst)
#
# logging.info("best num of feats: {}".format(n_feats[n_feats.index(best_cut_off)]))
#
# best_X = full_X[n_feats.index(best_cut_off)]
#
# logging.info("doing final SHAP")
# X_to_tes t= shap_feat_select(X_to_test, _dir=rad_file_path.parent, n_feat_cutoff=n_feats[n_feats.index(best_cut_off)],
#                             org=org)
#
# logging.info("Splitting data into 80-20 split")
# x_train, x_test, y_train, y_test = train_test_split(X_to_test, X_to_test.index, test_size=0.2, random_state=0)
#
#
# fig, ax = plt.subplots(figsize=[50, 50])
# sns.pairplot(x_train)
#
# ax.figure.axes[-1].yaxis.label.set_size(22)
#
# plt.tight_layout()
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
#
# plt.savefig(str(rad_file_path.parent) + "/pair_plot.png")
# plt.close()
#
# for x in range(20):
#     model = CatBoostClassifier(iterations=best_ntrees, random_state=x)
#     model.fit(x_train, y_train)
#     model_file_name = str(org_dir / "finalised_model.sav")
#     pickle.dump(model, open(model_file_name, 'wb'))
#     # logging.info("final model score: {}". format(model.score(x_test, y_test)))
#     logging.info("final ROC AUC score: {}".format(accuracy_score(y_test, model.predict(x_test))))
#     # logging.info("decision path: {}".format(model.decision_path(x_test)))
#     # logging.info("leaf indices: {}".format(model.apply(x_test)))