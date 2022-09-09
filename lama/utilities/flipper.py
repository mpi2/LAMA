import SimpleITK as sitk
import os
from pathlib import Path
from lama import common

def main(target_dir: Path = os.getcwd()):
    target_dir = Path("E:/try_emap_to_SD/permuted")
    volpaths = common.get_file_paths(target_dir)
    print('flipping')

    for path in volpaths:
        print(path)
        loader = common.LoadImage(path)
        vol = loader.img
        flipped_vol = sitk.Flip(vol, [True, True, True])
        #pa = sitk.PermuteAxesImageFilter()
        #pa.SetOrder([1,0,2])
        #flipped_vol = pa.Execute(vol)


        flipped_vol.SetDirection([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])
        flipped_vol.SetOrigin([0,0,0])
        sitk.WriteImage(flipped_vol, str(target_dir) + "/flipped_" + str(os.path.basename(path)), True)


if __name__ == '__main__':
    main()