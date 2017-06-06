#!/usr/bin/env python


import numpy as np
import SimpleITK as sitk
import sys
import os
from os.path import abspath
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
            d = abspath(d)
            print('\nGetting volumes from: ', d)
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
        z, y, x = zyx
        try:
            for k in self.data:
                out[k] = self.data[k][z, y, x]
            print(out.values())
            # print '\n\n'
            # print [x for x in out.values()]
        except IndexError:
            print("\nIndex out of bounds\n")
        return [x for x in out.values()]

    def get_roi(self, zzyyxx):
        out = OrderedDict()
        (z1, z2), (y1, y2), (x1, x2) = zzyyxx
        print z1, z2, y1, y2, x1, x2
        try:
            for k in self.data:
                # print self.data[k]
                # print type(self.data[k])
                # print self.data[k].shape
                out[k] = np.mean(self.data[k][z1:z2, y1:y2, x1:x2])
                print(k, ": " + str(out[k]))
            print('\n\n')
            print([x for x in out.values()])
        except IndexError:
            print("\nIndex out of bounds\n")
        return [x for x in out.values()]


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("get pixel values of a bunch of images")
    parser.add_argument('-d', '--dirs', dest='dirs', nargs='+', help='A series of dirs with images', required=True)
    parser.add_argument('-r', '--roi', dest='roi', help='get means of an ROI', default=False, action='store_true')
    args = parser.parse_args()
    d = DataGetter(args.dirs)

    while True:
        if not args.roi:
            zyx_str = raw_input("Enter 'z y x' coordinates\n")
            zyx = [int(x) for x in zyx_str.split()]
            try:
                d.get_pixels(zyx)
            except IndexError:
                print("Index out of range\n\n")
        else:
            roi_str = raw_input("Enter 'z-z y-y x-x' roi coordinates\n")
            z1, z2, y1, y2, x1, x2 = [ int(x) for x in roi_str.split()]
            try:
                d.get_roi(((z1,z2),(y1,y2),(x1,x2)))
            except IndexError:
                print("Index out of range\n\n")


