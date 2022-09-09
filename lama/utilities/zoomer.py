import SimpleITK as sitk
import os
from pathlib import Path
from lama import common
import nrrd
import numpy as np
from scipy import ndimage

target_dir = Path("E:/try_emap_to_SD/zoom_z_axis")

volpaths = common.get_file_paths(target_dir)

scaled = "scaled"
print(target_dir)
print('zooming')

for path in volpaths:
    vol, v_head = nrrd.read(path)
    print(path)
    loader = common.LoadImage(path)
    img = loader.img

    zoomed = ndimage.zoom(vol, zoom=[1,1,1.6], mode='nearest', order=0)
    #

    nrrd.write(str(target_dir) + "/scaled_" + str(os.path.basename(path)), zoomed, header=v_head)
