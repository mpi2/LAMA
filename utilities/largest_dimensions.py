#!/usr/bin/env python

import sys
from os.path import join, dirname
sys.path.insert(0, join(dirname(__file__), '..'))
import common
import numpy as np

def run(folder):
    max_zyx = [0, 0, 0]
    for file_ in common.get_file_paths(folder):
        shape = common.img_path_to_array(file_).shape
        max_zyx = np.amax((max_zyx, shape), axis=0)

    print max_zyx


if __name__ == '__main__':

    folder = sys.argv[1]
    run(folder)