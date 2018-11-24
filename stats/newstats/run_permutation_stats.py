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


import numpy as np
import os
import common
import pandas as pd
from pathlib import Path
import sys
from logzero import logger as logging
sys.path.insert(0, Path(__file__).absolute() / '..')
import distributions
import p_thresholds


def get_all_files(root_dir, endswith):
    for root, dirs, files in os.walk(root_dir):
        for name in files:
            if name.endswith(endswith):
                yield os.path.join(root, name)


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


def get_organ_volume_data(root_dir: Path) -> Path:
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

            df = pd.read_csv(organ_vol_file, index_col=0)
            df['line'] = line_dir.name
            dataframes.append(df)

    # Write the concatenated organ vol file to single csv
    all_organs = pd.concat(dataframes)
    all_organs.index.name = 'vol'
    outpath = output_dir / common.ORGAN_VOLUME_CSV_FILE
    all_organs.to_csv(outpath)

    return all_organs


def get_staging_data(root_dir: Path) -> Path:
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

            df = pd.read_csv(staging_info, index_col=0)
            df['line'] = line_dir.name
            dataframes.append(df)

    # Write the concatenated staging info to the
    all_staging = pd.concat(dataframes)
    outpath = output_dir / common.STAGING_INFO_FILENAME
    all_staging.to_csv(outpath)

    return all_staging


def annotate_lines():
    """
    Iterate over all the mutant registration folder and add a csv of organ calls
    Returns
    -------

    """
    pass

def generate_mutant_distributions():
    """
    Generate mutant p-value distributions by regressiong organ volume on genotype + staging metric

    Returns
    -------

    """


def prepare_data(wt_organ_vol: pd.DataFrame,
                 wt_staging: pd.DataFrame,
                 mut_organ_vol: pd.DataFrame,
                 mut_staging: pd.DataFrame,
                 log_staging: bool=False,
                 log_dependent: bool=False) -> pd.DataFrame:
    """
    Do some preprocessing on the input DataFrames and concatenate into one

    Returns
    -------
    Dataframe ...

    """
    wt_staging.rename(columns={'value': 'crl'}, inplace=True)
    mut_staging.rename(columns={'value': 'crl'}, inplace=True)
    wt_staging.index = wt_staging.index.astype(str)

    # merge the organ vol
    organ_vols = pd.concat([wt_organ_vol, mut_organ_vol])

    # Drop any organ columns that has only zero values. These are the gaps in the label map caused by merging labels
    organ_vols = organ_vols.loc[:, (organ_vols != 0).any(axis=0)]

    # For the statsmodels linear mode to work, column names cannot start with a digid. Prefix with 'x'
    organ_vols.columns = [f'x{x}' if x.isdigit() else x for x in organ_vols.columns]

    if log_dependent:
        logging.info('logging dependent variable')
        log_res = np.log(organ_vols.drop(['line'], axis=1))
        line = organ_vols[['line']]
        organ_vols = pd.concat([log_res, line], axis=1)

    # Merge the staging data
    staging = pd.concat([wt_staging, mut_staging])

    if log_staging:
        print('logging staging metric')
        log_res = np.log(staging.drop(['line'], axis=1))  # TODO: not finished
        staging = pd.concat([log_res, staging['line']], axis=1)

    # Merge staging to the organvolume dataframe
    data = organ_vols.merge(right=staging, left_index=True, right_index=True, suffixes=['', '_delete'])
    data = data.drop(['line_delete'], axis=1)

    return data


def run(wt_dir: Path, mut_dir: Path, out_dir: Path, num_perms: int, log_dependent: bool=False):
    """
    Run the premutation-based stats pipeline

    Parameters
    ----------
    wt_dir
        Root of the wild type registration output
        This should contain an 'inputs' folder that contains a single baseline folder that contains multiuple specimen folders
    mut_dir
        Root of the mutant registration output
        This should contain 'inputs' folder that contains multiple mutant lines folder, each containing one or more mutant specimen folders
    out_dir
        Where to store the intermediate results of the permutation testing
    num_perms
        number of permutations to do
    """
    # Collate all the staging and organ volume data into csvs

    wt_staging = get_staging_data(wt_dir)
    mut_staging = get_staging_data(mut_dir)

    wt_organ_vol = get_organ_volume_data(wt_dir)
    mut_organ_vol = get_organ_volume_data(mut_dir)

    data = prepare_data(wt_organ_vol,
                        wt_staging,
                        mut_organ_vol,
                        mut_staging,
                        log_dependent)

    out_dir.mkdir(exist_ok=True)

    # Get the null distributions
    line_null, specimen_null = distributions.null(data, num_perms)

    null_line_pvals_file = out_dir / 'null_line_dist_pvalues.csv'
    null_specimen_pvals_file = out_dir / 'null_specimen_dist_pvalues.csv'

    # Write the null distributions to file
    line_null.to_csv(null_line_pvals_file)
    specimen_null.to_csv(null_specimen_pvals_file)

    # Get the alternative distribution
    line_alt, spec_alt = distributions.alternative(data)

    line_alt_pvals_file = out_dir / 'alt_line_dist_pvalues.csv'
    spec_alt_pvals_file = out_dir / 'alt_specimen_dist_pvalues.csv'

    # Write the alternative distributions to file
    line_alt.to_csv(line_alt_pvals_file)
    spec_alt.to_csv(spec_alt_pvals_file)

    line_organ_thresholds = p_thresholds.get_thresholds(line_null, line_alt)
    specimen_organ_thresholds = p_thresholds.get_thresholds(specimen_null, spec_alt)

    line_thresholds_path = out_dir / 'line_organ_p_thresholds.csv'
    spec_thresholds_path = out_dir / 'specimen_organ_p_thresholds.csv'

    line_organ_thresholds.to_csv(line_thresholds_path)
    specimen_organ_thresholds.to_csv(spec_thresholds_path)





if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser("Permutation-based stats")
    parser.add_argument('-w', '--wt_dir', dest='wt_dir', help='wildtype registration directory', type=argparse.FileType('w'),
                        required=True)
    parser.add_argument('-m', '--mut_dir', dest='mut_dir', help='mutant registration directory', type=argparse.FileType('w'),
                        required=True)
    parser.add_argument('-o', '--out_dir', dest='out_dir', help='permutation results directory', type=argparse.FileType('w'),
                        required=True)
    parser.add_argument('-n', '--num_perm', dest='num_perm', help='number of permutations to do', type=np.int,
                        required=False, default=1000)

    args = parser.parse_args()

    run(args.wt_dir, args.mut_dir, args.out_dir, args.num_perm)