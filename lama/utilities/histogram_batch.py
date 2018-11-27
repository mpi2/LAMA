#!/usr/bin/env python

import os
from os.path import join, relpath
import shutil
import SimpleITK as sitk
from matplotlib import pyplot as plt
import numpy as np

CSS_FILE = 'style.css'
FILE_SUFFIX = '.png'


XLIM = (0, 256)


def get_plot(im_path, label, numbins=256, remove_zeros=False, log=True):
    itk_img = sitk.ReadImage(im_path)
    array1 = sitk.GetArrayFromImage(itk_img)

    if array1.dtype in ('uint8', 'int8'):
        max_ = 256
    elif array1.max() <= 300:  # Don't use 255. If the  8bit data has been scaled/registered it might end up as 16bit with some values slightly above 255
        max_ = array1.max()
    else:
        max_ = 65536

    if remove_zeros:
       array1 = array1[array1 > 0]

    # bins = np.linspace(array1.min(), array1.max(), num=binsize)
    bins = np.linspace(array1.min(), max_, num=numbins)
    hist1, bins = np.histogram(array1, bins=numbins)
    if log:
        hist1 = np.log(hist1)
    width = 1 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.legend(loc='upper center',
          fancybox=True, shadow=True)
    plt.xlim(array1.min(), max_)
    plt.bar(center, hist1, color='blue', align='center', alpha=0.4, width=1.0)
    plt.legend()
    return plt


def batch(dir_, outdir, bins=None, remove_zeros=False, log=False):
    return
    file_paths = get_file_paths(dir_)
    for path in file_paths:
        basename = os.path.splitext(os.path.basename(path))[0]
        outpath = os.path.join(outdir, basename + '.png')
        plt = get_plot(path, basename, bins, remove_zeros, log=log)
        plt.savefig(outpath)
        plt.close()
    make_html(outdir, outdir)

def print_list_vertically(my_list):
    for i in my_list:
        print((str(i)))

def make_html(in_dir, out_dir):

    html = '<html><link rel="stylesheet" type="text/css" href="{}" /><body>\n'.format(CSS_FILE)
    stage = os.path.basename(in_dir)
    html += '<div class="title"> Hitograms for {}</div>\n'.format(stage)

    for img in os.listdir(in_dir):
        if not img.endswith(FILE_SUFFIX):
            continue
        base = os.path.basename(img)
        html += '<div class="specimen">'

        img_path = join(in_dir, img)

        html += '<div class="chart"><img class="chart_img" src="{}"/></div>\n'.format(relpath(img_path, out_dir))
        html += '<div class="title clear">{}</div><br><br>'.format(base)
        # Close specimen
        html += '</div>'
    html += '</body></html>'

    outfile = join(out_dir, 'histograms.html')
    with open(outfile, 'w') as html_fh:
        html_fh.write(html)
        html_fh.write(get_css())


def get_css():
    css = """<style>
     body{font-family: Arial}
    .title{width: 100%; padding: 20px; background-color: lightblue; margin-bottom: 20px}
    .chart_img{width: 400px}
    .specimen{float: left; padding-right: 20px}
    .clear{clear: both}
    .title{font-size: 12px}
    </style>"""
    return css

def get_file_paths(folder, extension_tuple=('.nrrd', '.tiff', '.tif', '.nii', '.bmp', 'jpg', 'mnc', 'vtk'), pattern=None):
    """
    Test whether input is a folder or a file. If a file or list, return it.
    If a dir, return all images within that directory.
    Optionally test for a pattern to sarch for in the filenames
    """
    if not os.path.isdir(folder):
        if isinstance(folder, str):
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
    parser = argparse.ArgumentParser("Make histograms from a folder of volumes")
    parser.add_argument('-i', '--indir', dest='folder', help='', required=True)
    parser.add_argument('-o', '--out', dest='out', help='out dir', required=True)
    parser.add_argument('-l', '--log', dest='log', help='show log intensities', required=False, default=False, action='store_true')
    parser.add_argument('-b', '--bins', dest='bins', help='num bins', required=False, default=256, type=int)
    parser.add_argument('-z', '--rz', dest='rmzero', help='remove zeros', required=False, default=False, action='store_true')
    args = parser.parse_args()

    if args.folder:
        batch(args.folder, args.out, bins=args.bins, remove_zeros=args.rmzero, log=args.log)
