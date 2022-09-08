from pathlib import Path
import nrrd
import os
import numpy as np
from logzero import logger as logging


def main():
    # NOTE DO NOT USE AS YOU WILL CLIP SHIT
    img_path = Path(
        "Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210713_ark_target/210602_C3H_avg_n18.nrrd")
    label_path = Path(
        "Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210713_ark_target/labelled_clip.nrrd")

    img, img_h = nrrd.read(img_path)
    l_clip, c_head = nrrd.read(label_path)

    img[l_clip != 0] = np.random.choice([38, 39, 40])

    img[(img < 0) | (img > 255)] = np.random.choice([38, 39, 40])

    fixed_name = "fixed_avg.nrrd"
    logging.info(str(os.path.dirname(img_path) + "/" + fixed_name))
    logging.info(np.amax(img), np.amin(img))

    nrrd.write(str(os.path.dirname(img_path) + "/" + fixed_name), img)


if __name__ == '__main__':
    main()
