import SimpleITK as sitk
import os
from pathlib import Path
from lama import common

def main(target_dir: Path = os.getcwd()):
    target_dir = Path("E:/220204_BQ_dataset/221218_BQ_run/new_inv_labels")
    volpaths = common.get_file_paths(target_dir)
    print('flipping')

    for path in volpaths:
        print(path)
        loader = common.LoadImage(path)
        vol = loader.img
        flipped_vol = sitk.Flip(vol, [False, True, False])
        #pa = sitk.PermuteAxesImageFilter()
        #pa.SetOrder([1,0,2])
        #flipped_vol = pa.Execute(vol)
        print(flipped_vol.GetDirection())

        direction = flipped_vol.GetDirection()
        #print(type(direction[6]))

        #direction[6] = direction[6] * -1


        flipped_vol.SetDirection([-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0])

        print(flipped_vol.GetDirection())

        #flipped_vol.SetOrigin([0,0,0])
        sitk.WriteImage(flipped_vol, str(target_dir) + "/flipped/" + str(os.path.basename(path)), True)


if __name__ == '__main__':
    main()