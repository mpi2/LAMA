#!/usr/bin/env python3

"""
This module takes a list of directories each containing a mutant line or baseline input
It will try to make a lock on the input file and remove a line/specimen to process after which it will releae the lock
This is to enable multiple machines to process the data concurrently
"""
import sys
import shutil
from pathlib import Path # if you haven't already done so
file = Path(__file__).resolve()
parent, root = file.parent, file.parents[1]
sys.path.append(str(parent))
import socket
import pandas as pd
from filelock import FileLock, Timeout
import run_lama


TIMEOUT = 10


def lama_job_runner(job_file, config_path, root_directory):

    job_file = Path(job_file).absolute()
    config_path = Path(config_path)
    config_name = config_path.name

    lock = FileLock(f'{job_file}.lock', timeout=TIMEOUT)

    write_index = False

    while True:

        try:
            with lock.acquire():
                df_jobs = pd.read_csv(job_file)
                print('is', lock.is_locked)

                # add Status columns to list
                if 'status' not in df_jobs:
                    write_index = True # The first time we write the file add a numeric index
                    df_jobs['status'] = 'to_run'
                    df_jobs['host'] = None

                # Get an unfinished job
                jobs_to_do = df_jobs[df_jobs['status'] == 'to_run']

                if len(jobs_to_do) < 1:
                    print("No more jobs left on jobs list")
                    raise SystemExit

                indx = jobs_to_do.index[0]

                try:
                    dir_ = jobs_to_do.at[indx, 'dir']
                except KeyError:
                    print('p')

                df_jobs.at[indx, 'status'] = 'running'
                df_jobs.at[indx, 'host'] = socket.gethostname()
                df_jobs.to_csv(job_file, index=write_index)
                write_index = False

                dest_config_path = Path(root_directory) / dir_ / config_name
                # Move the config to the line/baseline input folder
                shutil.copy(config_path, dest_config_path)


        except Timeout:
            sys.exit('Timed out', socket.gethostname())

        try:
            print(f'trying {dir_}')
            run_lama.RegistrationPipeline(dest_config_path)

        except Exception as e:
            with lock.acquire():
                df_jobs = pd.read_csv(job_file, index_col=0)
                df_jobs.at[indx, 'status'] = 'failed'
                df_jobs.to_csv(job_file)

        else:
            with lock.acquire():
                df_jobs = pd.read_csv(job_file, index_col=0)
                df_jobs.at[indx, 'status'] = 'complete'
                df_jobs.to_csv(job_file)


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser("Schedule LAMA jobs")
    parser.add_argument('-j', '--job_list', dest='job_file', help='file_with jobs list watch for new jobs',
                        required=True)
    parser.add_argument('-c', '--config', dest='config', help='_pheno_detect.yaml config file',
                        required=True)
    parser.add_argument('-r', '--root_dir', dest='root_dir', help='The root directory containing the input folders',
                        required=True)


    args = parser.parse_args()
    lama_job_runner(Path(args.job_file), Path(args.config), Path(args.root_dir))


