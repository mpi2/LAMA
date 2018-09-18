#!/usr/bin/env python3

"""
This module submits jobs via qsub to the grid. The harwell grid as 16 nodes and each node will be used to run a baselines specimen or mutant line

work in progress

"""
import sys
import shutil
from pathlib import Path # if you haven't already done so
file = Path(__file__).resolve()
parent, root = file.parent, file.parents[1]
sys.path.append(str(parent))
import socket
import pandas as pd
from filelock import SoftFileLock, FileLock, Timeout
import run_lama
import os
import job_runner

TIMEOUT = 10


def run_on_grid():
    """
    Run LAMA on a grid using qsub to submit jobs.

    The LAMA docker must be uploaded to the local docker server first
    Returns
    -------

    """
    pass

def run_on_mutlipe_servers():
    """
    If LAMA is installed on multiple machines, ssh into each machine and set off the jobrunner
    Returns
    -------

    """
    pass


if __name__ == '__main__':

    # import argparse
    #
    # parser = argparse.ArgumentParser("Schedule LAMA jobs")
    # parser.add_argument('-j', '--job_list', dest='job_file', help='file_with jobs list watch for new jobs',
    #                     required=True)
    # parser.add_argument('-c', '--config', dest='config', help='_pheno_detect.yaml config file',
    #                     required=True)
    # parser.add_argument('-r', '--root_dir', dest='root_dir', help='The root directory containing the input folders',
    #                     required=True)

    run_on_grid(Path()


