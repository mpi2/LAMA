#!/usr/bin/env python

import SimpleITK as sitk
import numpy as np
import sys
import os
import csv


def img_stats(dir_, outfile=None):
    imgs = GetFilePaths(dir_)

    header = ['name', 'size', 'shape (z, y, x)', 'min', 'max', 'mean', 'spacing', 'origin', 'dtype', 'num 1s']

    rows = [header]

    for img_name in imgs:

        img_path = os.path.join(dir_, img_name)
        im = sitk.ReadImage(img_path)
        a = sitk.GetArrayFromImage(im)
        basename = os.path.basename(img_path)
        row = [basename, a.size, '"{}"'.format(a.shape), a.min(), a.max(), np.mean(a), '"{}"'.format(im.GetSpacing()),
               '"{}"'.format(im.GetOrigin()), a.dtype, a[a==1].size]
        rows.append(row)

    if outfile:
        with open(outfile, 'wb') as fh:
            writer = csv.writer(fh)
            writer.writerows(rows)
    else:
        for row in rows:
            print ', '.join([str(x) for x in row])

def GetFilePaths(folder, extension_tuple=('.nrrd', '.tiff', '.tif', '.nii', '.bmp', 'jpg', 'mnc', 'vtk', 'bin'), pattern=None):
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

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("Get stats on images in a folder")
    parser.add_argument('-i', '--indir', dest='dir', help='directory with images', required=True)
    parser.add_argument('-o', '--out_file', dest='out_file', help='Optional path to csv file', default=False)
    args = parser.parse_args()
    img_stats(args.dir, args.out_file)
