import nrrd
import random
import SimpleITK as sitk
from lama import common
import os
import pandas as pd
from pathlib import Path


def randomise_file_list(_dir):

    o_dir = Path(_dir.parent / "output")
    os.mkdir(o_dir)

    file_list = [file_name for file_name in common.get_file_paths(_dir)]
    print(file_list)
    # randomise list
    random.shuffle(file_list)

    file_df = pd.DataFrame({"name": file_list, "num": range(len(file_list))})

    print(file_df)

    file_df.to_csv(o_dir/"results.csv")
    for i, file_name in enumerate(file_list):
        i_name = o_dir / (str(i) + ".nrrd")
        os.rename(file_name,i_name)



def main():
    _dir = Path("//anufiles.anu.edu.au/anu/jcsmr/ArkellLab/Lab Members/Amrit/2022_BL6_cohort/BL6_full_no_genotype")
    randomise_file_list(_dir)

if __name__ == '__main__':
    main()