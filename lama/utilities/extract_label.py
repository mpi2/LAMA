from pathlib import Path
import nrrd
import os
from lama import common
import SimpleITK as sitk
from scipy import ndimage
import numpy as np

def main(target_dir, labs_of_interest: list =[17]):

    label_paths = [spec_path for spec_path in common.get_file_paths(target_dir) if ('inverted_labels' in str(spec_path))]
    #labs_of_interest = list(labs_of_interest)

    #labs_of_interest = [float(i) for i in labs_of_interest] if len(labs_of_interest) > 1 else labs_of_interest
    print(labs_of_interest)
    # Just get the label and write them
    for path in label_paths:
        label, l_head = nrrd.read(path)
        #print((~np.isin(label, labs_of_interest)))

        label[~np.isin(label, labs_of_interest)] = 0
        #label[label not in labs_of_interest] = 0
        #print(str(os.path.basename(path)))
        file_name = target_dir / "inverted_labels" / str(os.path.basename(path))
        nrrd.write(str(file_name), label)