#!/usr/bin/env python


from utilities.dev import read_minc
import os
import SimpleITK as sitk

def minc_to_nrrd(indir, outdir):
    paths = get_paths(indir)
    for path in paths:
        arr = read_minc.minc_to_numpy(path)
        im = sitk.GetImageFromArray(arr)
        basename = os.path.splitext(os.path.basename(path))[0]
        outpath = os.path.join(outdir, basename + '.nrrd')
        sitk.WriteImage(im, outpath)


def get_paths(indir):
    """
    Test whether input is a folder or a file. If a file or list, return it.
    If a dir, return all images within that directory.
    Optionally test for a pattern to sarch for in the filenames
    """
    if not os.path.isdir(indir):
        return False
    else:
        paths = []
        for root, _, files in os.walk(indir):
            for filename in files:
                if filename.lower().endswith('mnc'):
                    paths.append(os.path.join(root, filename))
        return paths




if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("batch convert minc to nrrd")
    parser.add_argument('-i', '--input', dest='input_dir', help='dir containing mincs',
                        required=True)
    parser.add_argument('-o', '--out', dest='out_dir', help='Dir to put nrrd into', required=True)

    args = parser.parse_args()
    minc_to_nrrd(args.input_dir, args.out_dir)