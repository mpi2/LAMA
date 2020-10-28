#!/usr/bin/env python3

import nrrd
from pathlib import Path
import itertools
import numpy as np
from lama import common


def convert_16_bit_to_8bit(indirs, outdir, clobber: bool):

    paths = [common.get_file_paths(Path(x)) for x in indirs]
    paths = itertools.chain(*paths)

    for inpath in paths:
        print(f'doing {inpath.name}')
        arr, header = nrrd.read(str(inpath))

        if arr.dtype not in (np.uint16, np.int16):
            print(("skipping {}. Not 16bit".format(inpath.name)))
            arr_cast = arr

        else:
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

        if clobber:
            outpath = inpath
        else:
            outpath = Path(outdir) / inpath.name
        nrrd.write(str(outpath), arr_cast, header=header)


def main():
    import argparse
    parser = argparse.ArgumentParser("Rescale 16 bit images to 8bit")
    parser.add_argument('-i', dest='indirs', help='dir with vols to convert. Will include subdirectories', required=True,  nargs='*')
    parser.add_argument('-o', dest='outdir', help='dir to put vols in. Omit to overwtrite source and use --clobber', required=False,
                        default=None)
    parser.add_argument('--clobber', dest='clobber', help='Overwrite inputs!', action='store_true',
                        default=None)

    args = parser.parse_args()

    if args.outdir and args.clobber:
        raise SystemExit("Use either --clobber OR -o")
    if args.outdir:
        outdir = Path(args.outdir)
    else:
        outdir = None

    # indirs = [Path(x) for x in args.indirs]

    convert_16_bit_to_8bit(args.indirs, outdir, args.clobber)

if __name__ == '__main__':
    main()

