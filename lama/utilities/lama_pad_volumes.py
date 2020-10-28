#! /usr/bin/env python3

"""
Pad one or more folder of 3D volumes to either
    the largest dimensions across all volumes
    or
    the to a set size

Examples
--------

# To pad a series of volumes inplace
$ lama_pad_volumes.py -i dir1 dir2 --clobber

# To pad a series of volumes and write to new directory
$ lama_pad_volumes.py -i dir1 dir2 --outdir new_dir

# To pad a series of volumes to a defined shape (xyz)
$ lama_pad_volumes.py -i dir1 dir2 --outdir new_dir -d 100 200 300


"""


from os.path import basename
import sys
from typing import Iterable, Tuple
from pathlib import Path

import SimpleITK as sitk
import nrrd
from logzero import logger as logging

from lama import common


def get_largest_dimensions(indirs: Iterable[Path]) -> Tuple[int]:
    max_dims = None

    for dir_ in indirs:

        volpaths = common.get_file_paths(dir_)

        for path in volpaths:
            im = sitk.ReadImage(str(path))
            dims = im.GetSize()

            if not max_dims:
                max_dims = dims
            else:
                max_dims = [max(d[0], d[1]) for d in zip(dims, max_dims)]

    return max_dims


def pad_volumes(indirs: Iterable[Path], max_dims: Tuple, outdir: Path, clobber: bool, filetype: str='nrrd'):
    """
    Pad volumes, masks, labels. Output files will have same name as original, but be in a new output folder

    Parameters
    ----------
    indirs
        one or more directories containing volumes to pad (Will search subdirectories for volumes)
    max_dims
        dimensions to pad to (z, y, x)
    outdir
        path to output dir
    """

    if clobber and outdir:
        print('Specifiy either --clobber or an output dir (-o)')
        return
    if not clobber and not outdir:
        print('Specifiy either --clobber or an output dir (-o)')
        return

    if not max_dims:
        max_dims = get_largest_dimensions(indirs)

    print(f'Zero padding to {max_dims}')

    outdir = outdir

    for dir_ in indirs:
        dir_ = Path(dir_)

        if clobber:
            result_dir = dir_
        else:
            result_dir = outdir / dir_.name
            result_dir.mkdir(exist_ok=True, parents=True)

        volpaths = common.get_file_paths(dir_)

        # print('Padding to {} - {} volumes/masks:'.format(str(max_dims), str(len(volpaths))))
        # pad_info = Dict()

        for path in volpaths:

            if clobber:
                outpath = path
            else:
                outpath = result_dir / path.name

            loader = common.LoadImage(path)
            vol = loader.img
            if not vol:
                logging.error('error loading image for padding: {}'.format(loader.error_msg))
                sys.exit()
            vol_dims = vol.GetSize()

            # The voxel differences between the vol dims and the max dims
            diffs = [m - v for m, v in zip(max_dims, vol_dims)]

            # How many pixels to add to the upper bounds of each dimension, divide by two and round down to nearest int
            upper_extend = [d // 2 for d in diffs]

            # In case of differnces that cannot be /2. Get the remainder to add to the lower bound
            remainders = [d % 2 for d in diffs]

            # Add the remainders to the upper bound extension to get the lower bound extension
            lower_extend = [u + r for u, r in zip(upper_extend, remainders)]

            # if any values are negative, stop. We need all volumes to be the same size
            for ex_val in zip(lower_extend, upper_extend):

                if ex_val[0] < 0 or ex_val[1] < 0:
                    msg = ("\ncan't pad images\n"
                           "{} is larger than the specified volume size\n"
                           "Current vol size:{},\n"
                           "Max vol size: {}"
                           "\nCheck the 'pad_dims' in the config file\n".format(basename(path), str(vol_dims),
                                                                                str(max_dims)))

                    logging.error(msg)
                    raise common.LamaDataException(msg)

            # Pad the volume. New pixels set to zero
            padded_vol = sitk.ConstantPad(vol, upper_extend, lower_extend, 0)
            padded_vol.SetOrigin((0, 0, 0))
            padded_vol.SetSpacing((1, 1, 1))


            sitk.WriteImage(padded_vol, str(outpath), True)
            # pad_info['data'][input_basename]['pad'] = [upper_extend, lower_extend]
    print('Finished padding')


def main():
    import argparse

    if len(sys.argv) < 2:
        print(__doc__)
        return
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--indirs', dest='indirs', help='directory with images', nargs='*', required=True)
    parser.add_argument('-o', '--outdir', dest='outdir', help='where to put padded images', default=None)
    parser.add_argument('-c', '--clobber', dest='clobber', help='force overwriting of input volumes', action='store_true', default=None)
    parser.add_argument('-d', '--new_dims', dest='new_dims', nargs=3, type=int, help='x y z to pad to (eg 100 150 300)',
                        required=False, default=False)

    args = parser.parse_args()

    if args.outdir:
        outdir = Path(args.outdir)
    else:
        outdir = None

    indirs = [Path(x) for x in args.indirs]
    pad_volumes(indirs, args.new_dims, outdir, args.clobber)


if __name__ == '__main__':
    main()

