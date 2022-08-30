# assumes CT image nrrds with minimum value of -1024
# folder with mask nrrd files
# folder with mask nrrd files
# assumes CT mask nrrds with 0's and 1's
# assumes the shape of each patient data (iamge and mask) are different - therefore,
# this will pad all images and masks to the size of the largest
# does not preform any interpolation to isotrpic voxels or any normalization
# only saves the image and mask, therefore the metadata and pixel spacing is lost

import nrrd # pip install pynrrd # probably better performance with sitk
import numpy as np
import h5py
from lama import common
from pathlib import Path
import random
import pandas as pd
import os

# input
# both folder should have the same number of files in the same order.. obviously..
# folder with image nrrd files

_dir =  Path('Z:/jcsmr/ROLab/Experimental data/Radiomics/Workflow design and trial results/Kyle Drover analysis/220827_pytorch-contouring/')
image_nrrd_folder = Path('Z:/jcsmr/ROLab/Experimental data/Radiomics/Workflow design and trial results/Kyle Drover analysis/220827_pytorch-contouring/imgs')

# folder with mask nrrd files
mask_nrrd_folder = Path('Z:/jcsmr/ROLab/Experimental data/Radiomics/Workflow design and trial results/Kyle Drover analysis/220827_pytorch-contouring/tumour_respaced')
# output
output = 'output'
# dataset name
dataset = "someName"


n_train = 50
n_validate = 12

# substitute value if larger
def addLargest(value,variable):
    if value>variable:
        variable=value
    return variable

# replace globs with common.getImages.

images = common.get_file_paths(image_nrrd_folder)
masks = common.get_file_paths(mask_nrrd_folder)


img_dirs = [x for x in images]
mask_dirs = [x for x in masks]

# sanity
assert(len(images)==len(masks))

dirs = pd.DataFrame({'imgs': img_dirs, 'labels': mask_dirs})

# shuffle rows in dataframe and reset index


dirs = dirs.sample(frac=1).reset_index(drop=True)


for i in range(n_train+n_validate):
    spec_dir = dirs.iloc[i]

    image_nrrd, img_h = nrrd.read(spec_dir['imgs'])

    mask_nrrd, mask_h = nrrd.read(spec_dir['labels'])

    if i < n_train:
        file5 = str(_dir) + '/train/' + f"{os.path.splitext(spec_dir['imgs'].name)[0]}.h5"
    elif i < (n_train + n_validate):
        file5 = str(_dir) + '/validate/' + f"{os.path.splitext(spec_dir['imgs'].name)[0]}.h5"
    else:
        break
    with h5py.File(file5, 'w') as f5:
        f5.create_dataset("raw", dtype=np.float32, data=image_nrrd)
        f5.create_dataset("label", dtype=np.uint8, data=mask_nrrd)
