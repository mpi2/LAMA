"""
The entry point to running permutation-based statistics.

Currently moving the permutation scripts from the work_scripts repo


How this module will work

Currently, this module works only with organ volume data. The voxel-based methods are currently too big to do this.
Think about including voxel-based stats in the future

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
import common
import struct
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, Path(__file__).absolute() / '..')
import permute
def generate_baseline_data():
    """

    Returns
    -------

    """

def prepare_mutant_data(mutant_dir, organ_info, out_path, stats_folder_name, specimen_level=False):
    """
    Permuter takes all the baseline from all specimens in csv files.
    This function just aggregates the mutant p value data from individual files so data so permuter can use it

    format of new dataframe should be specimen_id(index), label_num, label_num
    """
    organ_vol_pvalue_dfs = []

    def add_csv(path, first):
        if 'het' in path.lower():  # This removes lline
            print(path)
            return

        # extract p-values
        df = pd.read_csv(path, index_col=0)

        organ_vol_pvalue_dfs.append(df[['p']].T)  # .T to make organs as columns

        cols = df[['p']].T.columns
        if first:
            test = set(cols)
            first = False

    first = True

    for line_dir in [join(mutant_dir, x) for x in os.listdir(mutant_dir)]:
        if not os.path.isdir(line_dir): continue

        if not specimen_level:  # Get pvalues for the line level results
            organ_vol_csv = join(line_dir, 'output', 'stats', stats_folder_name, 'inverted_organ_volumes_LinearModel_FDR5%.csv')
            if not os.path.isfile(organ_vol_csv):
                print(f"Can't find {organ_vol_csv}")
                continue
            add_csv(organ_vol_csv, first)
        else:
            speciemn_level_dir = join(line_dir, 'output', 'stats', stats_folder_name, 'specimen_calls')
            for spec_csv in get_all_files(speciemn_level_dir, '.csv'):
                add_csv(spec_csv, first)

    # merge into one dataframe
    df_result = pd.concat(organ_vol_pvalue_dfs, axis=0)

    df_result = df_result.sort_index(axis=1)


    df_result.to_csv(out_path)

# Move to common
def get_all_files(root_dir, endswith):
    for root, dirs, files in os.walk(root_dir):
        for name in files:
            if name.endswith(endswith):
                yield os.path.join(root, name)


def get_organ_volume_data(root_dir: Path):
    """
    Given a root registration dorectory, collate all the organ volume csvs into one file
    Parameters
    ----------
    root_dir

    Returns
    -------

    """

    # extract p-values

    output_dir = root_dir / 'output'

    dataframes = []

    for line_dir, specimen_dir in iterate_over_specimens(output_dir):

            organ_vol_file = specimen_dir / 'output' / common.ORGAN_VOLUME_CSV_FILE

            if not organ_vol_file.is_file():
                raise FileNotFoundError(f'Cannot find organ volume file {organ_vol_file}')

            df = pd.read_csv(organ_vol_file)
            df['line'] = line_dir.name
            dataframes.append(df)
    # Write the concatenated organ vol file to signle csv
    all_staging = pd.concat(dataframes)
    all_staging.to_csv(output_dir / common.ORGAN_VOLUME_CSV_FILE)



def iterate_over_specimens(reg_out_dir: Path):

    if not reg_out_dir.is_dir():
        raise FileNotFoundError(f'Cannot find output directory {reg_out_dir}')

    for line_dir in reg_out_dir.iterdir():

        if not line_dir.is_dir():
            continue

        for specimen_dir in line_dir.iterdir():
            if not specimen_dir.is_dir():
                continue

            yield line_dir, specimen_dir


def get_staging_data(root_dir: Path):
    """
    Given a root directory for either baselines or mutants, collate all the staging data into a single csv

    Parameters
    ----------
    root_dir: Directory containing registration output.
        eg baselines

    """
    output_dir = root_dir / 'output'

    dataframes = []

    for line_dir, specimen_dir in iterate_over_specimens(output_dir):

            staging_info = specimen_dir / 'output' / common.STAGING_INFO_FILENAME

            if not staging_info.is_file():
                raise FileNotFoundError(f'Cannot find staging info file {staging_info}')

            df = pd.read_csv(staging_info)
            df['line'] = line_dir.name
            dataframes.append(df)
    # Write the concatenated staging info to the
    all_staging = pd.concat(dataframes)
    all_staging.to_csv(output_dir / common.STAGING_INFO_FILENAME)


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

    # get_staging_data(wt_dir)
    # get_staging_data(mut_dir)

    get_organ_volume_data(wt_dir)
    get_organ_volume_data(mut_dir)

    # null_distrbutions = permute.permuter()
    # mutant_per_organ_pvalues = permute.permuter()
    # per_organ_thresholds = get_thresholds(null_distrbutions, mutant_per_organ_pvalues)




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