from pathlib import Path
from lama import common
import SimpleITK as sitk
import numpy as np






def main():
    #path = Path("E:/220901_emap_e12_atlas/TS20_EMA76_reference.nrrd")

    #loader = common.LoadImage(path)
    #img = loader.img

    # get the otsu mask

    #inv = sitk.InvertIntensity(img)


    #sitk.WriteImage(inv, "E:/220901_emap_e12_atlas/TS20_EMA76_reference_inv.nrrd")

    path = Path("E:/220901_emap_e12_atlas/TS20_EMA76_reference_inv.nrrd")

    loader = common.LoadImage(path)
    img = loader.img

    euler_transform = sitk.Euler3DTransform()
    euler_transform.SetRotation(angleX=90)













if __name__ == '__main__':
    main()