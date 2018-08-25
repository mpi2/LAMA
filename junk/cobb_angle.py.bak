#!/usr/bin/env python

"""

"""

import numpy as np
import os
from os.path import join, basename, abspath

def cobb_angle(indir):
    dirs = [join(indir, x) for x in os.listdir(indir)]
    indir = abspath(indir)
    for dir_ in dirs:
        vtk_file = [join(indir, dir_, x) for x in os.listdir(dir_) if x.endswith('.vtk')][0]
        points = parse_vtk(vtk_file)
        thor1, thor2, lumb1, lumb2 = points
        CD = thor2 - thor1
        EF = lumb2 - lumb1
        angle_rad = angle_between(CD, EF)
        angles_degrees = np.degrees(angle_rad)
        print("{}\t{}".format(basename(dir_), angles_degrees))


def angle_between(v1, v2):
    """ Finds angle between two vectors `v1` and `v2`. """
    return np.arccos(np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))


def parse_vtk(infile):
    points = []
    with open(os.path.abspath(infile), 'r') as fh:
        for line in fh:
            try:
                float(line.split()[0])
            except ValueError:
                continue
            x, y, z = [float(x) for x in line.split()]
            points.append(np.array([x, z]))
    return points

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("icalculate Cobb angles from VTK files")
    parser.add_argument('-i', '--input', dest='input_dir', help='input_vtk', required=True)
    args = parser.parse_args()
    cobb_angle(args.input_dir)



