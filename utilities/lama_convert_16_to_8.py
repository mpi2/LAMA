#!/usr/bin/env python3

import SimpleITK as sitk
from pathlib import Path
import numpy as np
from lama import common


def convert_16_bit_to_8bit(indir, outdir):

    clobber = True if not outdir else False

    paths = common.get_file_paths(Path(indir))

    for inpath in paths:
        img = sitk.ReadImage(str(inpath))
        arr = sitk.GetArrayFromImage(img)

        if arr.dtype not in (np.uint16, np.int16):
            print(("skipping {}. Not 16bit".format(inpath.name)))
            continue

        if arr.max() <= 255:
            print(("16bit image but with 8 bit intensity range {}".format(inpath.name)))
            arr_cast = arr.astype(np.uint8)
            
        # Fix the negative values, which can be caused by  the registration process. therwise we end up with hihglights
        # where there should be black
        else:
            if arr.dtype == np.int16:
                # transform to unsigned range
                print('unsigned short')
                negative_range = np.power(2, 16) / 2
                arr += negative_range
            # Do the cast
            arr2 = arr/256
            arr_cast = arr2.astype(np.uint8)
            print((arr_cast.min(), arr_cast.max()))

        out_img = sitk.GetImageFromArray(arr_cast)

        if clobber:
            outpath = inpath
        else:
            outpath = Path(outdir) / inpath.name
        sitk.WriteImage(out_img, str(outpath), True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser("Rescale 16 bit images to 8bit")
    parser.add_argument('-i', dest='indir', help='dir with vols to convert. Will include subdirectories', required=True)
    parser.add_argument('-o', dest='outdir', help='dir to put vols in. omit to overwtrite source', required=False, default=None)
    args = parser.parse_args()
    convert_16_bit_to_8bit(args.indir, args.outdir)

