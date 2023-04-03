import numpy as np
from lama.utilities.prep_for_man_valid import resample, rotate
from lama import common
import nrrd
import SimpleITK as sitk
from pathlib import Path
import os
from logzero import logger as logging

def main():
    # NOTE DO NOT USE AS YOU WILL CLIP SHIT
    img_path = Path("E:/try_emap_to_SD/test_pad/220909_sd_rev_reg/target/TS20_EMA76_reference.nrrd")
    img = common.LoadImage(img_path).img
    img_list = []


    #padded_img = sitk.ConstantPad(img, pad_amount, pad_amount,np.min(sitk.GetArrayFromImage(img)).astype(float))
    #img_list.append(sitk.GetArrayFromImage(padded_img))

    img_list.append(sitk.GetArrayFromImage(img))

    rot_name = "rot2_TS20_EMA76_reference.nrrd"
    logging.info(str(os.path.dirname(img_path)+"/"+rot_name))
    rot_avg = rotate(img_list,0, -0.05, 0)
    #rot_avg_v2 = rotate(rot_avg, 0, 0, np.deg2rad(-180))

    nrrd.write(str(os.path.dirname(img_path)+"/"+rot_name), rot_avg[0])

if __name__ == '__main__':
    main()