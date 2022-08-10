from pathlib import Path
from lama import common
import SimpleITK as sitk
import os


def main():
    _dir = Path("E:/220809_tcp_paper_dataset/4D_atlas_nrrds/")
    masked = 'masked'
    for i, path in enumerate(common.get_file_paths(_dir)):
        loader = common.LoadImage(path)
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
        dilate.SetKernelRadius([0.9, 0.9, 0.9])
        dilate.SetKernelType(sitk.sitkBall)
        mask = dilate.Execute(mask)

        sitk.WriteImage(mask, str( _dir.parent / masked / os.path.basename(path)))


if __name__ == '__main__':
    main()