from lama.utilities.prep_for_man_valid import resample, rotate
import nrrd
import SimpleITK as sitk
from pathlib import Path
import os
from logzero import logger as logging

def main():
    # NOTE DO NOT USE AS YOU WILL CLIP SHIT
    img_path = Path("Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210713_ark_target/deformable_8.nrrd")
    img, img_h = nrrd.read(img_path)
    img_list = []
    img_list.append(img)

    rot_name = "rot_atlas.nrrd"
    logging.info(str(os.path.dirname(img_path)+"/"+rot_name))
    rot_avg = rotate(img_list, -0.075, 0.025, 0.025)

    nrrd.write(str(os.path.dirname(img_path)+"/"+rot_name), rot_avg[0])

if __name__ == '__main__':
    main()