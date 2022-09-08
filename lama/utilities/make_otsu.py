from pathlib import Path
from lama import common
import SimpleITK as sitk
import os
from body_radiomics_normaliser import fill_image
from logzero import logger as logging
import numpy as np

def main():
    _dir = Path("E:/try_emap_to_SD/TS20_EMA76_reference_inv_orient.nrrd")
    print(_dir)
    masked = 'masked'
    #for i, path in enumerate(common.get_file_paths(_dir)):


    loader = common.LoadImage(_dir)
    img = loader.img

    # get the otsu mask
    Otsu = sitk.OtsuThresholdImageFilter()

    inv_mask = Otsu.Execute(img)
    mask = sitk.InvertIntensity(inv_mask, 1)

    mask = sitk.ConnectedComponent(mask != mask[0, 0, 0])

    # sitk.WriteImage(seg, os.path.join(output, name + "_all_connected.nrrd"))
    mask = sitk.RelabelComponent(mask)
    mask = mask == 1
    # sitk.WriteImage(seg, os.path.join(output, name + "_largest_connected.nrrd"))

    # lets see if dilate with a tight kernal fixes getting stupid dots everywhere.
    dilate = sitk.BinaryDilateImageFilter()
    dilate.SetKernelRadius([1, 1, 1])
    dilate.SetKernelType(sitk.sitkBall)
    mask = dilate.Execute(mask)
    #npa  = sitk.GetArrayFromImage(mask)
    #logging.info("fill holes in first orientation")
    #npa_hole_filled = fill_image(npa)

    #logging.info("fill holes in second orientation")
    #npa_hole_filled = fill_image(npa_hole_filled, roll=1)

    #logging.info("fill holes in third orientation")
    #npa_hole_filled = fill_image(npa, roll=0)

    #transposed = np.transpose(npa_hole_filled, axes=(0, 2, 1))

    # Turn np array to image
    filler = sitk.VotingBinaryIterativeHoleFillingImageFilter()
    filler.SetMaximumNumberOfIterations(1000)
    filled = filler.Execute(mask)
    filler.CopyInformation(mask)


    sitk.WriteImage(filled, str( _dir.parent / masked / os.path.basename(_dir)))


if __name__ == '__main__':
    main()