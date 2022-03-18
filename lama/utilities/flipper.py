import SimpleITK as sitk
import os
from pathlib import Path
from lama import common

def main(target_dir: Path = os.getcwd()):

    volpaths = common.get_file_paths(target_dir)
    print('flipping')

    for path in volpaths:
        loader = common.LoadImage(path)
        vol = loader.img
        flipped_vol = sitk.Flip(vol, [True, False, False])

        flipped_vol.SetDirection([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])
        flipped_vol.SetOrigin([0,0,0])
        sitk.WriteImage(flipped_vol, str(path), True)