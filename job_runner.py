#!/usr/bin/env python3

"""
This module takes a list of directories each containing a mutant line or baseline input
It will try to make a lock on the input file and remove a line/specimen to process after which it will releae the lock
This is to enable multiple machines to process the data concurrently
"""
import sys
from pathlib import Path
file = Path(__file__).resolve()
parent, root = file.parent, file.parents[1]
sys.path.append(str(parent))

import shutil
from stats import run_lama_stats
import socket
import pandas as pd
from filelock import FileLock, Timeout
import run_lama
from datetime import datetime
from logzero import logger as logging

TIMEOUT = 10

allowed_types = ['registration',
                 'stats']

JOBFILE_NAME = 'lama_jobs.csv'

def prepare_inputs(jobs_file: Path, root_dir: Path, is_mutants:bool=False):
    """
    The inputs will be in a seperate folder. This function splits them out into individual subdirectories
    for lama to work on.

    It also copies the config file across

    It creates a joblist csv file for use with lama_job_runner

    Parameters
    ----------
    jobfile_path: Path to save job file to
    root_dir: the root project directory
    is_mutants: if True search the folder for individual line sub folders

    """

    output_dir = root_dir / 'output'
    output_dir.mkdir(exist_ok=True)

    jobs_entries = []

    if not is_mutants:
        # Look for inputs folder where all the baseline volumes should be
        inputs_dir = root_dir / 'inputs'
        if not inputs_dir.is_dir():
            raise FileNotFoundError("Cannot find baseline 'inputs' directory")

        # Create output folders for each volume and copy over volume
        for vol in inputs_dir.iterdir():
            vol_name = vol.stem
            specimen_inputs_dir = output_dir / vol_name / 'inputs'  # Lama will for volumes in 'inputs'
            specimen_inputs_dir.mkdir(exist_ok=True, parents=True)
            shutil.copy(vol, specimen_inputs_dir)

            # Create a job entry. Dir will be the specimen directory relative to the jobs file
            rel_path_to_specimen_root_dir = specimen_inputs_dir.parent.relative_to(jobs_file.parent)
            jobs_entries.append([rel_path_to_specimen_root_dir, 'to_run', '_', '_', '_'])

        jobs_df = pd.DataFrame.from_records(jobs_entries, columns=['dir', 'status', 'host', 'start_time', 'end_time'])

        jobs_df.to_csv(jobs_file)


def lama_job_runner(config_path: Path,
                    root_directory: Path,
                    type_: str='registration',
                    is_mutants: bool=False):

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
    is_mutants: If true, we are using a set of mutants, which will be split into line subfolders

    Returns
    -------

    """

    if type_ not in allowed_types:
        raise ValueError(f'{type_} must be one of {str(allowed_types)}')

    if not config_path.is_file():
        raise  FileNotFoundError(f"can't find config file {config_path}")

    job_file = root_directory / JOBFILE_NAME

    prepare_inputs(job_file, root_directory, is_mutants)

    config_name = config_path.name

    lock = FileLock(f'{job_file}.lock', timeout=TIMEOUT)

    while True:

        try:
            with lock:

                df_jobs = pd.read_csv(job_file, index_col=0)

                # Get an unfinished job
                jobs_to_do = df_jobs[df_jobs['status'] == 'to_run']

                if len(jobs_to_do) < 1:
                    print("No more jobs left on jobs list")
                    raise SystemExit(0)

                indx = jobs_to_do.index[0]

                dir_ = jobs_to_do.at[indx, 'dir']

                df_jobs.at[indx, 'status'] = 'running'

                df_jobs.at[indx, 'start_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

                df_jobs.at[indx, 'host'] = socket.gethostname()

                df_jobs.to_csv(job_file)

                if type_ == 'registration':
                    # Copy the config into each project directory
                    dest_config_path = Path(root_directory) / dir_ / config_name

                elif type_ == 'stats':
                    # copy the folder into output/stats/
                    # May not exist so make
                    dest_dir = Path(root_directory) / dir_ / 'output' / 'stats'
                    dest_dir.mkdir(exist_ok=True, parents=True)
                    dest_config_path = dest_dir / config_name

                # elif type_ == 'glcms':
                #     dest_dir = Path(root_directory) / dir_ / 'output' / 'glcm'
                #     dest_dir.mkdir(exist_ok=True, parents=True)
                #     dest_config_path = dest_dir / config_name


                shutil.copy(config_path, dest_config_path)

        except Timeout:
            sys.exit('Timed out' + socket.gethostname())

        try:
            print(f'trying {dir_}')

            if type_ == 'registration':

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

    parser = argparse.ArgumentParser("Schedule LAMA jobs")

    parser.add_argument('-c', '--config', dest='config', help='lama.yaml config file',
                        required=True)
    parser.add_argument('-r', '--root_dir', dest='root_dir', help='The root directory containing the input folders',
                        required=True)
    parser.add_argument('-t', '--type', dest='type_', help=f'one of {str(allowed_types)}', choices=allowed_types,
                        required=True)
    parser.add_argument('-m', '--mutants', dest='mutants', help='If -m is used we have a folder of mutants in subfolders',
                        required=False, type=bool, default='store_false')

    args = parser.parse_args()

    lama_job_runner(Path(args.job_file), Path(args.config), Path(args.root_dir), args.type_, args.mutants)


