from pathlib import Path
import nrrd
import os
from lama import common
import SimpleITK as sitk
from scipy import ndimage

def main(target_dir, labs_of_interest: list =[17]):

    label_paths = [spec_path for spec_path in common.get_file_paths(target_dir) if ('inverted_labels' in str(spec_path))]

    # Just get the label and write them
    for path in label_paths:
        label, l_head = nrrd.read(path)
        label[label not in labs_of_interest] = 0
        print(str(os.path.basename(path)))
        file_name = target_dir / "inverted_labels" / str(os.path.basename(path))
        nrrd.write(str(file_name), label)