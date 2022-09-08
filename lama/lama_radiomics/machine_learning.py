

from lama.lama_radiomics import feature_reduction
from lama import common

from filelock import SoftFileLock, Timeout
import socket
from datetime import datetime
import sys
import signal

from sklearn.ensemble import RandomForestClassifier, StackingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import  MLPClassifier
from sklearn.metrics import accuracy_score, auc, roc_auc_score, f1_score, matthews_corrcoef, precision_score, recall_score

from logzero import logger as logging
from pathlib import Path
import pandas as pd
import os

#from mpi4py import MPI
import sys

from joblib import Parallel, delayed

JOBFILE_NAME = 'ml_jobs.csv'

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


def make_ml_jobs_file(jobs_file: Path, file_paths: list):
    """
    Creates a joblist csv file for use with the radiomics pipeline.
    Searches for all images paths and creates a job file

    Parameters
    ----------
    jobfile_path: Path to save job file to
    root_dir: the root project directory
    is_mutants: if True search the folder for individual line sub folders

    """
    # output_dir = root_dir / 'radiomics_output'
    # output_dir.mkdir(exist_ok=True)

    jobs_entries = []
    # get each file path
    for i, csv_path in enumerate(file_paths):

        rel_path_to_org_input = str(csv_path.relative_to(jobs_file.parent))
        jobs_entries.append([rel_path_to_org_input, 'to_run', '_', '_', '_'])

    jobs_df = pd.DataFrame.from_records(jobs_entries, columns=['job', 'status', 'host', 'start_time', 'end_time'])

    jobs_df.to_csv(jobs_file)
    return True



def ml_job_runner(org_dir):
    '''i
    Performs the pyradiomic calculations


    Parameters
    ----------
    target_dir



    Returns
    -------

    '''


    # get org csv files

    org_dir = Path(org_dir)
    names = common.get_file_paths(org_dir, extension_tuple=".csv")



    jobs_file_path = org_dir / JOBFILE_NAME
    lock_file = jobs_file_path.with_suffix('.lock')
    lock = SoftFileLock(lock_file)

    if not os.path.exists(jobs_file_path):
        logging.info("Creating a job-file for ml")
        make_ml_jobs_file(jobs_file_path, names)
        logging.info("Job_file_created")

    df_jobs = pd.read_csv(jobs_file_path, index_col=0)

    # execute parallelisation:
    while True:
        try:
            with lock.acquire(timeout=60):

                # Get an unfinished job
                jobs_to_do = df_jobs[df_jobs['status'] == 'to_run']
                if len(jobs_to_do) < 1:
                    logging.info("No more jobs left on jobs list")

                    # error trap for processes that hung
                    logging.info("checking for hung jobs")
                    fin_jobs = df_jobs[df_jobs['status'] == 'completed']
                    t_last_job_run = fin_jobs['start_time'].max()

                    # scan start time of running jobs - if they started before the latest
                    # completed job - it hung
                    hung_jobs = df_jobs[(df_jobs['status'] == 'running') & (df_jobs['start_time'] < t_last_job_run)]
                    if len(hung_jobs) > 0:
                        logging.info("Hung jobs found - rerunning")
                        jobs_to_do = hung_jobs
                    else:
                        break
                indx = jobs_to_do.index[0]

                org_csv_path = Path(org_dir) / (jobs_to_do.at[indx, 'job'])




                df_jobs.at[indx, 'status'] = 'running'
                df_jobs.at[indx, 'start_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                df_jobs.at[indx, 'host'] = socket.gethostname()

                df_jobs.to_csv(jobs_file_path)
        except Timeout:
            sys.exit('Timed out' + socket.gethostname())

        # try:
        logging.info(f'trying {org_csv_path}')
        # get the organ file and number
        org_df = pd.read_csv(org_csv_path)

        # didn't remove org before
        org = org_df['org'][0]

        # perform feature reduction on a single organ

        feature_reduction.main(org_df, org, org_dir)

        # except Exception as e:
        #    if e.__class__.__name__ == 'KeyboardInterrupt':
        #        logging.info('terminating')
        #        sys.exit('Exiting')

        #    status = 'failed'
        #    print(e)
        #    logging.exception(e)

        status = 'complete'

        with lock:
            df_jobs = pd.read_csv(jobs_file_path, index_col=0)
            df_jobs.at[indx, 'status'] = status
            df_jobs.at[indx, 'end_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            df_jobs.to_csv(jobs_file_path)

    logging.info('Exiting job_runner')
    return True


def main():
    import argparse

    parser = argparse.ArgumentParser("Run RF models for prediction")


    parser.add_argument('-i', '--input_file', dest='indirs', help='radiomics file', required=True,
                        type=str)

    args = parser.parse_args()



    rad_file_path = Path(args.indirs)
    X = pd.read_csv(str(rad_file_path))
    #run feature reduction in parallel
    Parallel(n_jobs=-1, verbose=2)(delayed(feature_reduction.main(X, org=org, rad_file_path=rad_file_path))(org) for org in X['org'].unique())



if __name__ == '__main__':
    main()