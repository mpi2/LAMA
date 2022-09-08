"""Just roughly normalises intensities of volumes and straightens each embryo to perform some sort of
 manual phenotyping"""
from lama.img_processing import normalise
from logzero import logger as logging
from lama import common
import os
import nrrd
from pathlib import Path
from scipy import ndimage
import numpy as np
import SimpleITK as sitk


def get_images(dir, s):
    img_list = []
    spec_name_list = []
    int_paths = common.get_file_paths(dir)

    # enumerating for speed only
    for i, img_path in enumerate(int_paths):
        img, img_h = nrrd.read(img_path)
        # only get heatmap vals inside of the mask + padding
        img = img[s[0].start:s[0].stop,
              s[1].start:s[1].stop,
              s[2].start:s[2].stop]
        spec_name_list.append(os.path.splitext(img_path.name)[0])

        img_list.append(img)
    return img_list, spec_name_list


def resample(image, transform):
    """
  This function resamples (updates) an image using a specified transform
  :param image: The sitk image we are trying to transform
  :param transform: An sitk transform (ex. resizing, rotation, etc.
  :return: The transformed sitk image
  """
    reference_image = image
    interpolator = sitk.sitkBSpline
    default_value = 0
    return sitk.Resample(image, reference_image, transform,
                         interpolator, default_value)


def rotate(vols, x, y, z):
    # TODO: fix so it doesn't clip
    logging.info(f"rotating vol with manual rotation {x, y, z}")
    rotated = []
    for vol in vols:
        print(np.shape(vol))
        fixed = sitk.GetImageFromArray(vol.astype(np.uint8), isVector=False)
        rigid_euler = sitk.Euler3DTransform()
        rigid_euler.SetRotation(x, y, z)
        rigid_euler.SetTranslation((0, 0, 0))
        # set center as midpoints
        rigid_euler.SetCenter([coord // 2 for coord in np.shape(vol)[::-1]])
        # rigid_euler.TransformPoint([point for point in grid for grid in img])
        mov = resample(fixed, rigid_euler)
        rotated.append(sitk.GetArrayFromImage(mov).astype(np.uint8))

    return rotated


def main():
    wt_dir = Path(
        "Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210521_vis_anal/wt")
    mut_dir = Path(
        "Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210521_vis_anal/non_wt")

    mask, mask_h = nrrd.read(
        "Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210423_g_by_e_stand_out/210415_g_by_e_anal/target/stats_mask.nrrd")

    pop_avg, pop_h = nrrd.read(
        "Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210423_g_by_e_stand_out/210415_g_by_e_anal/target/210224_pop_avg_deformable_8.nrrd")

    s = ndimage.find_objects(mask)[0]

    # get the images
    wt_imgs, wt_names = get_images(wt_dir, s)

    mut_imgs, mut_names = get_images(mut_dir, s)

    int_norm = normalise.IntensityMaskNormalise()

    # normalise the images
    int_norm.add_reference(wt_imgs)

    int_norm.normalise(mut_imgs)

    int_norm.normalise(wt_imgs)

    # manually orient the embryos
    wt_imgs = rotate(wt_imgs, -0.0635, 0.02, 0.01)

    mut_imgs = rotate(mut_imgs, -0.0635, 0.02, 0.01)

    # write files
    logging.info('writing files')
    for i, vol in enumerate(mut_imgs):
        nrrd.write(mut_names[i] + ".nrrd", vol)

    for i, vol in enumerate(wt_imgs):
        nrrd.write(wt_names[i] + ".nrrd", vol)


if __name__ == '__main__':
    main()
