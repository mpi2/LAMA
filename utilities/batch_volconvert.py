#!/usr/bin/env python
from scipy.stats.distributions import make_dict

import SimpleITK as sitk
import sys
import os
import argparse


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


def convert(indir, outdir, type_, cast, makedirs):

    vols = GetFilePaths(indir)

    if len(sys.argv) < 4:
        sys.exit('Use sitk to convert volume types. Usage eg: <in_dir> <out_dir> <file extension eg tif>')

    for volpath in vols:
        out = outdir
        if makedirs:
            dir_ = os.path.splitext(os.path.basename(volpath))[0]
            individual_dir = os.path.join(outdir, dir_)
            os.mkdir(individual_dir)
            out = individual_dir
        vol = sitk.ReadImage(volpath)
        if cast:
            if cast == 'uchar':
                c = sitk.sitkUInt8
            if cast == 'char':
                c = sitk.sitkInt8
            if cast == 'ushort':
                c = sitk.sitkUInt16
            if cast == 'short':
                c = sitk.sitkInt16
            if cast == 'float':
                c = sitk.sitkFloat32
            vol = sitk.Cast(vol, c)
        volid = os.path.splitext(os.path.basename(volpath))[0]
        newpath = os.path.join(out, volid + '.' + type_)
        sitk.WriteImage(vol, newpath)

if __name__ == '__main__':
    parser = argparse.ArgumentParser("The MRC Harwell phenotype detection tool")
    parser.add_argument('-in', dest='indir', help='dir with vols to convert', required=True)
    parser.add_argument('-out', dest='outdir', help='dir to put vols in', required=True)
    parser.add_argument('-t', dest='type', help='file extension for new vols', required=True)
    parser.add_argument('-c', dest='cast', help='char, short, float. If not selected, no casting takes palce',
                        default=None)
    parser.add_argument('-d', dest='makedirs', help='put each vol imto seperate folders', action='store_true')
    args = parser.parse_args()
    convert(args.indir, args.outdir, args.type, args.cast, args.makedirs)
