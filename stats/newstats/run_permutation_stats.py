"""
The entry point to running permutation-based statistics.

Currently moving the permutation scripts from the work_scripts repo


How this module will work

Currently, this module works only with organ volume data. The voxel-based methods are currently too big to do this.
u think about including voxel-based stats in the future

generate_baseline_data
generate_mutant_data
generate_mutant_staging_data
-----------------------------------------------
These functions search the registration output folders for the csvs that contain the organ volumes and collate into
single csvs for baseline and mutant. Also
"""


import subprocess as sub
import numpy as np
import os
import struct
import pandas as pd
from pathlib import Path
from . import permute

def generate_baseline_data():
    """

    Returns
    -------

    """

def generate_mutant_data() -> Path:
    pass

def generate_mutant_staging_data() -> Path:
    pass

def generate_baseline_staging_data() -> Path:
    pass

def make_null_distribution() -> Path:
    pass

def make_mutant_distribution() -> Path:
    pass

def get_thresholds() -> pd.DataFrame:
    pass

def annotate_lines():
    """
    Iterate over all the mutant registration folder and add a csv of organ calls
    Returns
    -------

    """
    pass


def run(wt_dir: Path, mut_dir: Path, num_perms: int):
    # Call all the functions to get the data
    null_distrbutions = permute.permuter()
    mutant_per_organ_pvalues = permute.permuter()
    per_organ_thresholds = get_thresholds(null_distrbutions, mutant_per_organ_pvalues)




if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser("Permutation-based stats")
    parser.add_argument('-w', '--wt_dir', dest='wt_dir', help='mutant registration directory', type=argparse.FileType('w'),
                        required=True)
    parser.add_argument('-m', '--mut_dir', dest='mut_dir', help='mutant registration directory', type=argparse.FileType('w'),
                        required=True)
    parser.add_argument('-n', '--num_perm', dest='num_perm', help='number of permutations to do', type=np.int,
                        required=False, default=1000)
    args = parser.parse_args()
    run(args.wt_dir, args.mut_dir, args.num_perm)