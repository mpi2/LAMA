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

from filelock import FileLock, Timeout
from logzero import logger as logging
import pandas as pd
import toml

from lama.registration_pipeline import run_lama
from lama.registration_pipeline.validate_config import LamaConfigError

TIMEOUT = 10


JOBFILE_NAME = 'lama_jobs.csv'


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

    input_root_dir = root_dir / 'inputs'  # This will contain line folders or a baseline folder

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
                    root_directory: Path):

    """

    Parameters
    ----------
    config_path:
        path to registration config file:
    root_directory
        path to root directory. The folder names from job_file.dir will be appending to this path to resolve projject directories

    """
    if not config_path.is_file():
        raise FileNotFoundError(f"can't find config file {config_path}")

    root_directory = root_directory.resolve()

    job_file = root_directory / JOBFILE_NAME

    lock = FileLock(f'{job_file}.lock', timeout=TIMEOUT)

    # TODO: What happens when we run a second jobrunner but the first is still preapring inputs

    # If there's no job file, then this is the first instance of job_runner running
    # Preapre the data
    if not job_file.is_file():
        prepare_inputs(job_file, root_directory)

    config_name = config_path.name

    while True:

        try:
            with lock:

                df_jobs = pd.read_csv(job_file, index_col=0)

                # Get an unfinished job
                jobs_to_do = df_jobs[df_jobs['status'] == 'to_run']

                if len(jobs_to_do) < 1:
                    print("No more jobs left on jobs list")
                    raise SystemExit()

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
            print(f'trying {dir_}')
            run_lama.run(dest_config_path)

        except LamaConfigError as lce:
            status = 'config_error'
            logging.exception(f'There is a problem with the config\n{lce}')

        except Exception as e:
            if e.__class__.__name__ == 'KeyboardInterrupt':
                logging.info('terminating')
                sys.exit('Exiting')
            with lock.acquire():
                status = 'failed'
                logging.exception(e)

        else:
            with lock.acquire():
                status = 'complete'
        finally:
            df_jobs = pd.read_csv(job_file, index_col=0)
            df_jobs.at[indx, 'status'] = status
            df_jobs.at[indx, 'end_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            df_jobs.to_csv(job_file)


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser("Schedule LAMA jobs")

    parser.add_argument('-c', '--config', dest='config', help='lama.yaml config file',
                        required=True)
    parser.add_argument('-r', '--root_dir', dest='root_dir', help='The root directory containing the input folders',
                        required=True)
    args = parser.parse_args()

    lama_job_runner(Path(args.config), Path(args.root_dir))


