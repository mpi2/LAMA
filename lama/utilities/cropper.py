import SimpleITK as sitk
import os
from pathlib import Path
from lama import common
import nrrd
import numpy as np
from scipy import ndimage
from logzero import logger as logging

def main():


    target_dir = Path("E:/try_emap_to_SD/final")
    volpaths = common.get_file_paths(target_dir)

    cropped = "cropped"
    os.makedirs(str(target_dir/cropped), exist_ok=True)
    logging.info("Cropping")

    for path in volpaths:
        logging.info(f"Doing {os.path.basename(path)}")
        vol, v_head = nrrd.read(path)

        loader = common.LoadImage(path)
        img = loader.img

        # get the otsu mask
        Otsu = sitk.OtsuThresholdImageFilter()

        inv_mask = Otsu.Execute(img)
        mask = sitk.InvertIntensity(inv_mask, 1)

        mask = sitk.ConnectedComponent(mask != mask[0, 0, 0])

        #sitk.WriteImage(seg, os.path.join(output, name + "_all_connected.nrrd"))
        mask = sitk.RelabelComponent(mask)
        mask = mask == 1
        #sitk.WriteImage(seg, os.path.join(output, name + "_largest_connected.nrrd"))

        #lets see if dilate with a tight kernal fixes getting stupid dots everywhere.
        dilate = sitk.BinaryDilateImageFilter()
        dilate.SetKernelRadius([1,1,1])
        dilate.SetKernelType(sitk.sitkBall)
        mask = dilate.Execute(mask)

        #sitk.WriteImage(mask, str(target_dir / masked / os.path.basename(path)))

        mask_arr = sitk.GetArrayFromImage(mask)

        # get the bounding box of the mask

        s = ndimage.find_objects(mask_arr)[0]

        # Add some tight padding

        p = 3

        crop_vol = vol[s[2].start - p: s[2].stop + p,
                   s[1].start - p: s[1].stop + p,
                   s[0].start - p: s[0].stop + p]
        #
        #l_clip, c_head = nrrd.read(target_dir / clip / os.path.basename(path))

        #crop_vol[l_clip != 0] = np.random.choice([38,39,40])

        file_path = target_dir / cropped

        nrrd.write(str(file_path/ os.path.basename(path)), crop_vol, header=v_head)

if __name__ == '__main__':
    main()