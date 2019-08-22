#!/usr/bin/env python3

"""
This module takes a list of directories each containing a mutant line or baseline input
It will try to make a lock on the input file and remove a line/specimen to process after which it will releae the lock
This is to enable multiple machines to process the data concurrently
"""
import sys
import os
from pathlib import Path


# Bodge until I get imports working in Docker
lama_docker_dir = Path('/lama')
if lama_docker_dir.is_dir():
    print('setting lama path bodge')
    par = Path(__file__).parents[1].resolve()
    sys.path.append(str(par))
    print(sys.path)

import shutil
import socket
from datetime import datetime
import time

from filelock import SoftFileLock, Timeout
from logzero import logger as logging
import pandas as pd
import toml

from inspect import currentframe, getframeinfo

from lama.registration_pipeline import run_lama
from lama.registration_pipeline.validate_config import LamaConfigError



JOBFILE_NAME = 'lama_jobs.csv'


def linenum():
    cf = currentframe()
    return cf.f_back.f_lineno

def process_specimen(vol: Path, output_dir: Path, jobs_file: Path, jobs_entries):
    vol_name = vol.stem
    specimen_inputs_dir = output_dir / vol_name / 'inputs'  # Lama will look for volumes in 'inputs'
    specimen_inputs_dir.mkdir(exist_ok=True, parents=True)
    shutil.copy(vol, specimen_inputs_dir)

    # Create a job entry. Dir will be the specimen directory relative to the jobs file
    rel_path_to_specimen_root_dir = specimen_inputs_dir.parent.relative_to(jobs_file.parent)
    jobs_entries.append([rel_path_to_specimen_root_dir, 'to_run', '_', '_', '_'])


def prepare_inputs(jobs_file: Path, root_dir: Path):
    """
    The inputs will be in a seperate folder.
    This function splits them out into individual subdirectories based on line for lama to work on.

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

    logging.info('Copying input data')

    jobs_entries = []

    input_root_dir = root_dir / 'inputs'  # This will contain one or more line folders or a single baseline folder

    # Get the line subdirectories
    for line in input_root_dir.iterdir():
        if not line.is_dir():
            continue
        for vol in line.iterdir():
            line_outdir = output_dir / line.name
            process_specimen(vol, line_outdir, jobs_file, jobs_entries)

    jobs_df = pd.DataFrame.from_records(jobs_entries, columns=['dir', 'status', 'host', 'start_time', 'end_time'])

    jobs_df.to_csv(jobs_file)


def lama_job_runner(config_path: Path,
                    root_directory: Path,
                    is_first_instance=False):

    """

    Parameters
    ----------
    config_path:
        path to registration config file:
    root_directory
        path to root directory. The folder names from job_file.dir will be appending to this path to resolve project directories

    Notes
    -----
    This function uses a SoftFileLock for locking the job_file csv to prevent multiple instances of this code from
    processing the same line or specimen. A SoftFileLock works by creating a lock file, and the presence of this file
    prevents other instances from accessing it. We don't use FileLock (atlhough this is better) as it's not
    supported on nfs file systems.

    If this script terminates unexpectedly while it has a lock on the file, it will not be released and the file
    remains. Therefore before running this script, ensure no previous lock file is hanging around.

    If running multiple instances of lama_job_runner, pass in the argument --first_instance. This will instruct the
    current instance to setup the folders and delete any pr



    """

    if not config_path.is_file():
        raise FileNotFoundError(f"can't find config file {config_path}")

    root_directory = root_directory.resolve()

    job_file = root_directory / JOBFILE_NAME
    lock_file = job_file.with_suffix('.lock')
    lock = SoftFileLock(lock_file)
    init_file = root_directory / 'init'

    HN = socket.gethostname()

    if is_first_instance:
        if job_file.is_file():
            os.remove(job_file)
            # Delete any lockfile and startfile that might be present. As we are using SoftFileLock, it could be present if a previous
            # job_runner instance failed unexpectedly
        if lock_file.is_file():
            os.remove(lock_file)

        try:
            with lock.acquire(timeout=1):
                print(f'debug {HN}, {linenum()}')
                prepare_inputs(job_file, root_directory)
                init_file.touch() # This indicates to other instances that everything has been initialised
        except Timeout:
            print(f"Make sure lock file: {lock_file} is not present on running first instance")
            sys.exit()

    else:
        while True:

            if init_file.is_file():
                break
            else:
                time.sleep(5)
                print('Waiting for init file')

    config_name = config_path.name

    while True:

        try:
            print(f'debug {HN}, {linenum()}')

            with lock.acquire(timeout=60):

                print(f'debug {HN}, {linenum()}')
                # Create a lock then read jobs and add status to job file to ensure job is run once only.
                df_jobs = pd.read_csv(job_file, index_col=0)

                # Get an unfinished job
                jobs_to_do = df_jobs[df_jobs['status'] == 'to_run']

                if len(jobs_to_do) < 1:
                    print("No more jobs left on jobs list")
                    break

                indx = jobs_to_do.index[0]

                dir_ = jobs_to_do.at[indx, 'dir']

                df_jobs.at[indx, 'status'] = 'running'

                df_jobs.at[indx, 'start_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

                df_jobs.at[indx, 'host'] = socket.gethostname()

                df_jobs.to_csv(job_file)

                # Copy the config into each project directory
                dest_config_path = root_directory / dir_ / config_name

                if dest_config_path.is_file():
                    os.remove(dest_config_path)

                shutil.copy(config_path, dest_config_path)

                # rename the target_folder now we've moved the config
                c = toml.load(dest_config_path)

                target_folder = config_path.parent / c.get('target_folder')
                # Can't seem to get this to work with pathlib
                target_folder_relpath = os.path.relpath(target_folder, str(dest_config_path.parent))
                c['target_folder'] = target_folder_relpath

                with open(dest_config_path, 'w') as fh:
                    fh.write(toml.dumps(c))

        except Timeout:
            sys.exit('Timed out' + socket.gethostname())

        try:
            print(f'debug {HN}, {linenum()}')
            print(f'trying {dir_}')
            run_lama.run(dest_config_path)

        except LamaConfigError as lce:
            status = 'config_error'
            logging.exception(f'There is a problem with the config\n{lce}')
            sys.exit()

        except Exception as e:
            if e.__class__.__name__ == 'KeyboardInterrupt':
                logging.info('terminating')
                sys.exit('Exiting')

            status = 'failed'
            logging.exception(e)

        else:
            status = 'complete'

        finally:
            with lock:
                df_jobs = pd.read_csv(job_file, index_col=0)
                df_jobs.at[indx, 'status'] = status
                df_jobs.at[indx, 'end_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                df_jobs.to_csv(job_file)
    print('Exiting job_runner')
    return True


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser("Schedule LAMA jobs")

    parser.add_argument('-c', '--config', dest='config', help='lama.yaml config file',
                        required=True)
    parser.add_argument('-r', '--root_dir', dest='root_dir', help='The root directory containing the input folders',
                        required=True)
    parser.add_argument('-f', '--first_instance', dest='is_first_instance', help='Is this the first instance of job runner',
                    default=False)
    args = parser.parse_args()

    try:
        lama_job_runner(Path(args.config), Path(args.root_dir), args.is_first_instance)
    except pd.errors.EmptyDataError as e:
        logging.exception(f'poandas read failure {e}')


