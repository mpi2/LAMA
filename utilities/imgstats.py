#!/usr/bin/env python

import SimpleITK as sitk
import numpy as np
import sys
import os
import csv


def img_stats(dir_, outfile=None):
    imgs = os.listdir(dir_)

    header = ['name', 'size', 'shape', 'min', 'max', 'mean', 'spacing', 'origin', 'dtype']

    rows = [header]

    for img_name in imgs:

        img_path = os.path.join(dir_, img_name)
        im = sitk.ReadImage(img_path)
        a = sitk.GetArrayFromImage(im)
        basename = os.path.basename(img_path)
        row = [basename, a.size, '"{}"'.format(a.shape), a.min(), a.max(), np.mean(a), '"{}"'.format(im.GetSpacing()),
               '"{}"'.format(im.GetOrigin()), a.dtype]
        rows.append(row)

    if outfile:
        with open(outfile, 'wb') as fh:
            writer = csv.writer(fh)
            writer.writerows(rows)
    else:
        for row in rows:
            print ', '.join([str(x) for x in row])



if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("Get stats on images in a fodler")
    parser.add_argument('-d', '--dir', dest='dir', help='directory with images', required=True)
    parser.add_argument('-o', '--out_file', dest='out_file', help='Optional path to csv file', default=False)
    args = parser.parse_args()
    img_stats(args.dir, args.out_file)