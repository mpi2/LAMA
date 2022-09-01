

from lama.radiomics import feature_reduction


from sklearn.ensemble import RandomForestClassifier, StackingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import  MLPClassifier
from sklearn.metrics import accuracy_score, auc, roc_auc_score, f1_score, matthews_corrcoef, precision_score, recall_score


from pathlib import Path
import pandas as pd


from joblib import Parallel, delayed

def establish_model(X, stack: bool=False):
    # do feature selection:
    X = feat_select(X, 'accuracy')

    # set up training and test_data
    x_train, x_test, y_train, y_test = train_test_split(X, X.index, stratify=y, test_size=0.2, random_state=0)

    #train models


    rf = RandomForestClassifier()


    #for i, mod in enumerate(estimators):
    #    mod.fit(x_train, y_train)
    if stack: #LAMA radiomics
        knn = KNeighborsClassifier(n_neighbors=5, n_jobs=-1)
        mlp_adam = MLPClassifier(solver='adam')
        mlp_lbfgs = MLPClassifier(solver='lbfgs')
        estimators = [('knn', knn),
                      ('mlp_adam', mlp_adam),
                      ('mlp_lbfgs', mlp_lbfgs),
                      ('rf', rf)]

        stack_model = StackingClassifier(estimators=estimators, final_estimator=LogisticRegression, n_jobs=-1)
    else: #BQ Radiomics
        stack_model = rf

    stack_model.fit(x_train, y_train)
    y_train_predict = stack_model.predict(x_train)
    y_test_predict = stack_model.predict(x_test)


    metrics = [accuracy_score,
               auc,
               roc_auc_score,
               f1_score,
               matthews_corrcoef,
               precision_score,
               recall_score]
    for i, met in enumerate(metrics):
        train_acc = met(y_train, y_train_predict)
        test_acc = met(y_test, y_test_predict)
        print(met, train_acc, test_acc)


    return stack_model


def main():
    import argparse


    parser.add_argument('-i', '--input_file', dest='indirs', help='Raw NRRD directory', required=True,
                        type=str)

    args = parser.parse_args()



    rad_file_path = Path(args.indirs)
    X = pd.read_csv(str(rad_file_path))
    #run feature reduction in parallel
    Parallel(n_jobs=-1)(delayed(feature_reduction.main(X, org=org, rad_file_path=rad_file_path))(org) for org in X['org'].unique())


if __name__ == '__main__':
    main()