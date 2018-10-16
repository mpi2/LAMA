#!/usr/bin/env python3

"""
This module takes a list of directories each containing a mutant line or baseline input
It will try to make a lock on the input file and remove a line/specimen to process after which it will releae the lock
This is to enable multiple machines to process the data concurrently
"""
import sys
import shutil
from pathlib import Path
file = Path(__file__).resolve()
parent, root = file.parent, file.parents[1]
sys.path.append(str(parent))
from stats import run_lama_stats
import socket
import pandas as pd
from filelock import FileLock, Timeout
import run_lama
from datetime import datetime
from logzero import logger as logging

TIMEOUT = 10


def lama_job_runner(job_file: str, config_path: str, root_directory: str, type_: str='lama_registration'):
    """

    Parameters
    ----------
    job_file
        path to csv containing the job list
        columns:
            dir: name of folder with daat to run
    config_path:
        path to config file:
            either registration config or stats config
    root_directory
        path to root directory. The folder names from job_file.dir will be appending to this path to resolve projject directories
    type_:
        lama_regitration: run lama registration pipeline
        stats: run the stats pipeline

    Returns
    -------

    """

    allowed_types = ['lama_registration',
                     'stats']

    if type_ not in allowed_types:
        raise ValueError(f'{type_} must be one of {str(allowed_types)}')

    job_file = Path(job_file)
    if not job_file.is_file():
        raise FileNotFoundError(f"Can't find job file {job_file}")

    config_path = Path(config_path)
    if not config_path.is_file():
        raise  FileNotFoundError(f"can't find config file {config_path}")

    config_name = config_path.name

    lock = FileLock(f'{job_file}.lock', timeout=TIMEOUT)

    write_index = False

    while True:

        try:
            with lock:
                df_jobs = pd.read_csv(job_file)
                print('is', lock.is_locked)

                # add Status columns to list
                if 'status' not in df_jobs:
                    write_index = True  # The first time we write the file add a numeric index
                    df_jobs['status'] = 'to_run'
                    df_jobs['host'] = "_"
                    df_jobs['start_time'] = "_"
                    df_jobs['end_time'] = "_"

                # Get an unfinished job
                jobs_to_do = df_jobs[df_jobs['status'] == 'to_run']

                if len(jobs_to_do) < 1:
                    print("No more jobs left on jobs list")
                    raise SystemExit(0)

                indx = jobs_to_do.index[0]

                try:
                    dir_ = jobs_to_do.at[indx, 'dir']
                except KeyError:
                    print('index error')
                    raise SystemExit

                df_jobs.at[indx, 'status'] = 'running'
                df_jobs.at[indx, 'start_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                # df_jobs.at[indx, 'end_time'] = ""

                df_jobs.at[indx, 'host'] = socket.gethostname()  # TODO if the host column has nothing in it it is inititalised as numeric and dies here
                df_jobs.to_csv(job_file, index=write_index)
                write_index = False

                if type_ == 'lama_registration':
                    # Copy the config into each prject directory
                    dest_config_path = Path(root_directory) / dir_ / config_name

                elif type_ == 'stats':
                    # copy the folder into output/stats/
                    # May not exist so make
                    dest_dir = Path(root_directory) / dir_ / 'output' / 'stats'
                    dest_dir.mkdir(exist_ok=True, parents=True)
                    dest_config_path = dest_dir / config_name

                shutil.copy(config_path, dest_config_path)

        except Timeout:
            sys.exit('Timed out' + socket.gethostname())

        try:
            print(f'trying {dir_}')

            if type_ == 'lama_registration':

                run_lama.RegistrationPipeline(dest_config_path)

            elif type_ == 'stats':
                run_lama_stats.run(dest_config_path)

        except Exception as e:
            with lock.acquire():
                df_jobs = pd.read_csv(job_file, index_col=0)
                df_jobs.at[indx, 'status'] = 'failed'
                df_jobs.at[indx, 'end_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                df_jobs.to_csv(job_file)
                logging.exception(e)

        else:
            with lock.acquire():
                df_jobs = pd.read_csv(job_file, index_col=0)
                df_jobs.at[indx, 'status'] = 'complete'
                df_jobs.at[indx, 'end_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                df_jobs.to_csv(job_file)


if __name__ == '__main__':

    import argparse

    allowed_types = ['lama_registration',
                     'stats']

    parser = argparse.ArgumentParser("Schedule LAMA jobs")
    parser.add_argument('-j', '--job_list', dest='job_file', help='file_with jobs list watch for new jobs',
                        required=True)
    parser.add_argument('-c', '--config', dest='config', help='lama.yaml config file',
                        required=True)
    parser.add_argument('-r', '--root_dir', dest='root_dir', help='The root directory containing the input folders',
                        required=True)
    parser.add_argument('-t', '--type', dest='type_', help=f'one of {str(allowed_types)}',
                        required=True)

    args = parser.parse_args()
    lama_job_runner(Path(args.job_file), Path(args.config), Path(args.root_dir), args.type_)


