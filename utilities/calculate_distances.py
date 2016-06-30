import os
import re
import numpy as np
from scipy.spatial import distance
"""
Calculate distance between a series of points inverted using elastix
"""

in_dir = '/home/neil/sig/LAMA_results/E15.5/CBX4/090516_CBX4/mutant_run/output/organ_length_calcs/out/affine'
ext = '.vtk'


f1_f2 = []
f3_f4 = []
for path, subdirs, files in os.walk(in_dir):
    for name in files:
        if not name.endswith(ext):
            continue
        points_file = os.path.join(path, name)
        with open(points_file, 'rb') as fh:
            all_points = []
            for line in fh:
                try:
                    float(line[0])
                except:
                    continue
                points_str = line.strip().split(' ')
                points = np.array([float(i) for i in points_str])
                all_points.append(points)
        # f1_f2.append(np.linalg.norm(all_points[0] - all_points[1]))
        # f3_f4.append(np.linalg.norm(all_points[2] - all_points[3]))
        f1_f2.append(np.linalg.norm(distance.euclidean(all_points[0], all_points[1])))
        #f3_f4.append(np.linalg.norm(distance.euclidean(all_points[2], all_points[3])))

print 'f1-f2:\n {}'.format(f1_f2)
# print 'f3-f4:\n {}'.format(f3_f4)

