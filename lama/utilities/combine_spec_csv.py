from pathlib import Path
import nrrd
import os
import numpy as np
from logzero import logger as logging
from lama import common
from typing import Union, List, Tuple, Dict
import pandas as pd


def main():
    #replace Patth with your root registration results (ie. wildtype_and_mutant_data)
    root_dir = Path(
        "Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210801_g_by_e_anal/g_by_e_data")

    # big chunk of code creates a pandas dataset containing all the organ volumes for each
    # registration directory stored in the .csv file
    
    full_staging_data = pd.concat(
        [pd.read_csv(spec) for spec in common.get_file_paths(folder=root_dir, extension_tuple=".csv")
         if (common.STAGING_INFO_FILENAME in str(spec))],
        ignore_index=True) #Replace STAGING_INFO_FILENAME with ORGAN_INFO_FILENAME to get organ vols

    output = root_dir / "full_staging.csv"
    full_staging_data.to_csv(output)

    full_organ_data = pd.concat(
        [pd.read_csv(spec) for spec in common.get_file_paths(folder=root_dir, extension_tuple=".csv")
         if (common.ORGAN_VOLUME_CSV_FILE in str(spec))],
        ignore_index=True)  # Replace STAGING_INFO_FILENAME with ORGAN_INFO_FILENAME to get organ vols

    output = root_dir / "full_organs.csv"
    full_organ_data.to_csv(output)



if __name__ == '__main__':

    main()
