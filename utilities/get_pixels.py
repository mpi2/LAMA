#!/usr/bin/env python


import numpy as np
import SimpleITK as sitk
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import common
import tempfile
from collections import OrderedDict


class DataGetter(object):
    def __init__(self, dirs):
        self.dirs = dirs
        self.data = self.memorymap_data(dirs)

    def memorymap_data(self, dirs):
        imgs = OrderedDict()
        for d in dirs:
            print '\nGetting volumes from: ', d
            for imgpath in common.GetFilePaths(d):
                basename = os.path.basename(imgpath)
                arr = common.img_path_to_array(imgpath)
                t = tempfile.TemporaryFile()
                m = np.memmap(t, dtype=arr.dtype, mode='w+', shape=arr.shape)
                m[:] = arr
                imgs[basename] = m
        return imgs

    def get_pixels(self, zyx):
        out = OrderedDict()
        for k in self.data:
            out[k] = self.data[k][zyx]
            print k , ": " +  str(out[k])
        print '\n\n'
        print [x for x in out.values()]



if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("get pixel values of a bunch of images")
    parser.add_argument('-d', '--dirs', dest='dirs', nargs='+', help='A series of dirs with images', required=True)
    args = parser.parse_args()
    d = DataGetter(args.dirs)

    while True:
        zyx_str = raw_input("Enter 'z y x' coordinates\n")
        zyx = zyx_str.split()
        try:
            d.get_pixels(zyx)
        except IndexError:
            print "Index out of range\n\n"


