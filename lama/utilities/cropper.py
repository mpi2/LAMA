import SimpleITK as sitk
import os
from pathlib import Path
from lama import common
import nrrd
import numpy as np
from scipy import ndimage

target_dir = Path("Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210713_ark_target")

volpaths = common.get_file_paths(target_dir)
clip = "clip"
cropped = "cropped"
masked = "masked"
array = "array"
print('cropping')

for path in volpaths:
    good_vol, v_head = nrrd.read(path)

    print(np.amax(good_vol))
    if np.amax(good_vol) < 2:
        break
    # loader = common.LoadImage(path)
    # img = loader.img

    # get the otsu mask
    # Otsu = sitk.OtsuThresholdImageFilter()

    # inv_mask = Otsu.Execute(img)
    # mask = sitk.InvertIntensity(inv_mask, 1)

    # sitk.WriteImage(mask, str(target_dir / masked / os.path.basename(path)))

    mask_arr, m_head = nrrd.read(str(target_dir / masked / os.path.basename(path)))
    # mask_arr = sitk.GetArrayFromImage(mask)

    nrrd.write(str(target_dir / array / os.path.basename(path)), mask_arr, header=v_head)

    # get the bounding box of the mask

    s = ndimage.find_objects(mask_arr)[0]
    
    # Add some tight padding

    p = 3

    crop_vol = good_vol[s[0].start - p: s[0].stop + p,
               s[1].start - p: s[1].stop + p,
               s[2].start - p: s[2].stop + p]
    #
    l_clip, c_head = nrrd.read(target_dir / clip / os.path.basename(path))

    crop_vol[l_clip != 0] = np.random.choice([38,39,40])

    nrrd.write(str(target_dir / cropped / os.path.basename(path)), crop_vol, header=v_head)


