#!/usr/bin/env python3

"""
This module takes a directory containing one or more subdirectories each containing a mutant line or baseline inputs
It will try to make a lock on the input file and remove a line/specimen to process after which it will release the lock
This is to enable multiple machines to process the data concurrently.

"""
import sys
import os
from pathlib import Path
import logzero


# Bodge until I get imports working in Docker
lama_docker_dir = Path('/lama')
if lama_docker_dir.is_dir():
    par = Path(__file__).parents[1].resolve()
    sys.path.append(str(par))
    print(sys.path)

import shutil
import socket
from datetime import datetime

from filelock import SoftFileLock, Timeout
from logzero import logger as logging
import pandas as pd
import toml

from inspect import currentframe

from lama.registration_pipeline import run_lama
from lama.registration_pipeline.validate_config import LamaConfigError
from lama.common import cfg_load


JOBFILE_NAME = 'lama_jobs.csv'


def linenum():
    cf = currentframe()
    return cf.f_back.f_lineno


def make_jobs_file(jobs_file: Path, root_dir: Path):
    """
    Creates a joblist csv file for use with lama_job_runner.
    Searches for all images paths in subdirectories of root_dir/inputs
    and ...

    Parameters
    ----------
    jobfile_path: Path to save job file to
    root_dir: the root project directory
    is_mutants: if True search the folder for individual line sub folders

    """
    output_dir = root_dir / 'output'
    output_dir.mkdir(exist_ok=True)

    jobs_entries = []

    input_root_dir = root_dir / 'inputs'  # This will contain one or more line folders or a single baseline folder

    # Get the line subdirectories
    for line in input_root_dir.iterdir():
        if not line.is_dir():
            continue
        for vol_path in line.iterdir():

            # Create a job entry. Dir will be the specimen directory relative to the jobs file
            rel_path_to_specimen_input = str(vol_path.relative_to(root_dir))
            jobs_entries.append([rel_path_to_specimen_input, 'to_run', '_', '_', '_'])

    jobs_df = pd.DataFrame.from_records(jobs_entries, columns=['job', 'status', 'host', 'start_time', 'end_time'])

    jobs_df.to_csv(jobs_file)
    return True


def lama_job_runner(config_path: Path,
                    root_directory: Path,
                    make_job_file: bool=False,
                    log_level=None):

    """

    Parameters
    ----------
    config_path:
        path to registration config file:
    root_directory
        path to root directory. The folder names from job_file.dir will be appending to this path to resolve project directories
    make_job_file
        if true, just make the job_file that other instances can consume

    Notes
    -----
    This function uses a SoftFileLock for locking the job_file csv to prevent multiple instances of this code from
    processing the same line or specimen. A SoftFileLock works by creating a lock file, and the presence of this file
    prevents other instances from accessing it. We don't use FileLock (atlhough this is more robust) as it's not
    supported on nfs file systems. The advantage of SoftFileLock is you can create a lock file manually if
    you want to edit a job file manually while job_runner is running (make sure to delete after editing).

    If this script terminates unexpectedly while it has a lock on the file, it will not be released and the file
    remains. Therefore before running this script, ensure no previous lock file is hanging around.
    """
    if log_level:
        logzero.loglevel(log_level)

    if not config_path.is_file():
        raise FileNotFoundError(f"can't find config file {config_path}")

    root_directory = root_directory.resolve()

    job_file = root_directory / JOBFILE_NAME
    lock_file = job_file.with_suffix('.lock')
    lock = SoftFileLock(lock_file)
    # init_file = root_directory / 'init'

    HN = socket.gethostname()

    if make_job_file:

        # Delete any lockfile and job_file that might be present from previous runs.
        if job_file.is_file():
            os.remove(job_file)

        if lock_file.is_file():
            os.remove(lock_file)

        try:
            with lock.acquire(timeout=1):
                logging.info('Making job list file')
                make_jobs_file(job_file, root_directory)
                logging.info('Job file created!. You can now run job_runner from multiple machines')
                return

        except Timeout:
            logging.error(f"Make sure lock file: {lock_file} is not present on running first instance")
            sys.exit()

    config_name = config_path.name

    while True:

        try:
            # Create a lock then read jobs and add status to job file to ensure job is run once only.
            with lock.acquire(timeout=60):

                df_jobs = pd.read_csv(job_file, index_col=0)

                # Get an unfinished job
                jobs_to_do = df_jobs[df_jobs['status'] == 'to_run']

                if len(jobs_to_do) < 1:
                    logging.info("No more jobs left on jobs list")
                    break

                indx = jobs_to_do.index[0]

                vol = root_directory / (jobs_to_do.at[indx, 'job'])

                df_jobs.at[indx, 'status'] = 'running'

                df_jobs.at[indx, 'start_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

                df_jobs.at[indx, 'host'] = socket.gethostname()

                df_jobs.to_csv(job_file)

                # Make a project dir drectory for specimen
                # vol.parent should be the line name
                # vol.stem is the specimen name minus the extension
                spec_root_dir = root_directory / 'output' / vol.parent.name / vol.stem
                spec_input_dir = spec_root_dir / 'inputs'
                spec_input_dir.mkdir(exist_ok=True, parents=True)
                spec_out_dir = spec_root_dir / 'output'
                spec_out_dir.mkdir(exist_ok=True, parents=True)
                shutil.copy(vol, spec_input_dir)

                # Copy the config into the project directory
                dest_config_path = spec_root_dir / config_name

                if dest_config_path.is_file():
                    os.remove(dest_config_path)

                shutil.copy(config_path, dest_config_path)

                # rename the target_folder now we've moved the config
                c = cfg_load(dest_config_path)

                target_folder = config_path.parent / c.get('target_folder')
                # Can't seem to get this to work with pathlib
                target_folder_relpath = os.path.relpath(target_folder, str(dest_config_path.parent))
                c['target_folder'] = target_folder_relpath

                with open(dest_config_path, 'w') as fh:
                    fh.write(toml.dumps(c))

        except Timeout:
            sys.exit('Timed out' + socket.gethostname())

        try:
            logging.info(f'trying {vol.name}')
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
    logging.info('Exiting job_runner')
    return True


def main():
    import argparse

    parser = argparse.ArgumentParser("Schedule LAMA jobs")

    parser.add_argument('-c', '--config', dest='config', help='lama.yaml config file',
                        required=True)
    parser.add_argument('-r', '--root_dir', dest='root_dir', help='The root directory containing the input folders',
                        required=True)
    parser.add_argument('-m', '--make_job_file', dest='make_job_file', help='Run with this option forst to crate a job file',
                    action='store_true', default=False)
    args = parser.parse_args()

    try:
        lama_job_runner(Path(args.config), Path(args.root_dir), args.make_job_file)
    except pd.errors.EmptyDataError as e:
        logging.exception(f'poandas read failure {e}')


if __name__ == '__main__':
    main()


