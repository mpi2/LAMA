from pathlib import Path
import nrrd
import os
import numpy as np
from logzero import logger as logging
from lama import common
from typing import Union, List, Tuple, Dict
import pandas as pd


def main(root_dir):
    # big chunk of code creates a pandas dataset containing all the staging / organ volumes for each
    # registration directory stored in the .csv file
    print(root_dir)
    full_staging_data = pd.concat(
        [pd.read_csv(spec) for spec in common.get_file_paths(folder=root_dir, extension_tuple=".csv")
         if (common.STAGING_INFO_FILENAME in str(spec))],
        ignore_index=True)  # Replace STAGING_INFO_FILENAME with ORGAN_INFO_FILENAME to get organ vols

    # remove duplicates
    print(full_staging_data)
    full_staging_data.drop(columns=["value"], inplace=True)
    full_staging_data.replace('', np.nan, inplace=True)

    full_staging_data.dropna(subset=["staging"], inplace=True)
    print(full_staging_data)

    output = root_dir / "full_staging.csv"
    full_staging_data.to_csv(output, index=False)

    full_organ_data = pd.concat(
        [pd.read_csv(spec) for spec in common.get_file_paths(folder=root_dir, extension_tuple=".csv")
         if (common.ORGAN_VOLUME_CSV_FILE in str(spec))],
        ignore_index=True)  # Replace STAGING_INFO_FILENAME with ORGAN_INFO_FILENAME to get organ vols

    output = root_dir / "full_organs.csv"
    full_organ_data.to_csv(output)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser("Collect organ and staging files")
    parser.add_argument('-i', dest='indirs', help='Root registration directory',
                        required=True)

    args = parser.parse_args()
    main(args.indirs)
