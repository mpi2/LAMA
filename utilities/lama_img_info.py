#!/usr/bin/env python3

import SimpleITK as sitk
import numpy as np
import os
import pandas as pd


def img_stats(dir_, outfile=None):

    if not dir_:
        # Use current directory
        dir_ = os.getcwd()

    imgs = get_volume_paths(dir_)

    header = ['name', 'size', 'z', 'y', 'x', 'min', 'max', 'mean', 'spacing', 'origin', 'dtype']

    rows = []

    for img_path in imgs:

        im = sitk.ReadImage(img_path)
        a = sitk.GetArrayFromImage(im)
        basename = os.path.basename(img_path)
        row = [basename,
               a.size,
               int(a.shape[0]),
               int(a.shape[1]),
               int(a.shape[2]),
               a.min(),
               a.max(),
               np.around(np.mean(a), 3),
               '{}'.format(im.GetSpacing()),
               '{}'.format(im.GetOrigin()),
               a.dtype]

        rows.append(row)

    df = pd.DataFrame.from_records(rows, columns=header)

    print(df.to_string())

    print(f'\nmax dimensions (z: {df.z.max()}, y: {df.y.max()}, x: {df.x.max()})')

    if outfile:
        df.to_csv(outfile, index=0)


def get_volume_paths(folder, extension_tuple=('.nrrd', '.tiff', '.tif', '.nii', '.bmp', 'jpg', 'mnc', 'vtk', 'bin'), pattern=None):
    """
    Test whether input is a folder or a file. If a file or list, return it.
    If a dir, return all images within that directory.
    Optionally test for a pattern to sarch for in the filenames
    """
    if not os.path.isdir(folder):
        return False
    else:
        paths = []
        for root, _, files in os.walk(folder):
            for filename in files:
                if filename.lower().endswith(extension_tuple):
                    if pattern:
                        if pattern and pattern not in filename:
                            continue
                    #paths.append(os.path.abspath(os.path.join(root, filename))) #Broken on shared drive
                    paths.append(os.path.join(root, filename))
        return paths


def main():
    import argparse
    parser = argparse.ArgumentParser("Get stats on images in a folder")
    parser.add_argument('-i', '--indir', dest='dir', help='directory with images', default=None)
    parser.add_argument('-o', '--out_file', dest='out_file', help='Optional path to csv file', default=False)
    args = parser.parse_args()
    img_stats(args.dir, args.out_file)

if __name__ == '__main__':
    main()

