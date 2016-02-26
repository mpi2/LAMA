#!/usr/bin/env python

import os
import shutil
import SimpleITK as sitk
from matplotlib import pyplot as plt
import numpy as np




def get_plot(im_path, label, binsize, outpath, remove_zeros=False, log=False):
    itk_img = sitk.ReadImage(im_path)
    array1 = sitk.GetArrayFromImage(itk_img)
    hist1, bins = np.histogram(array1, bins=range(binsize))
    if remove_zeros:
        hist1 = np.delete(hist1, 0)#remove first histogram value (zeroes)
        bins = np.delete(bins, 0)
    if log:
        hist1 = np.log(hist1)
    width = 1* (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.legend(loc='upper center',
          fancybox=True, shadow=True)
    plt.bar(center, hist1, width=width, color='blue', align='center', alpha=0.4, linewidth=0)
    plt.legend()
    return plt


def single(img_path, bins, remove_zeros):

def batch(dir_, outdir, bins, remove_zeros):
    file_paths = get_file_paths(dir_)
    for path in file_paths:
        basename = os.path.splitext(os.path.basename(path))[0]
        outpath = os.path.join(outdir, basename + '.png')
        plt = get_plot(path, basename, bins, outpath, remove_zeros)
        plt.savefig(outpath)
        plt.close()

def print_list_vertically(my_list):
    for i in my_list:
        print str(i)

def get_file_paths(folder, extension_tuple=('.nrrd', '.tiff', '.tif', '.nii', '.bmp', 'jpg', 'mnc', 'vtk'), pattern=None):
    """
    Test whether input is a folder or a file. If a file or list, return it.
    If a dir, return all images within that directory.
    Optionally test for a pattern to sarch for in the filenames
    """
    if not os.path.isdir(folder):
        if isinstance(folder, basestring):
            return [folder]
        else:
            return folder
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


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser("Stats component of the phenotype detection pipeline")
    parser.add_argument('-d', '--dir', dest='folder', help='', default=False)
    parser.add_argument('-i', '--input', dest='image', help='', default=False)
    #parser.add_argument('-l', '--lab', dest='label', help='Label for pplot', required=True, type=str)
    parser.add_argument('-o', '--out', dest='out', help='out dir', required=False)
    parser.add_argument('-b', '--bins', dest='bins', help='Num bins', required=False, default=256, type=int)
    parser.add_argument('-z', '--rz', dest='rmzero', help='remove zeros', required=False, default=False, action='store_true')
    args = parser.parse_args()

    if args.folder
    batch(args.folder, args.out, args.bins, args.rmzero)
