#!/usr/bin/env python

"""
Convert a series of wt and mut volumes into a csv file that can be loaded into R for use in PhenStat
"""

import sys
from os.path import join
import os
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common
import numpy as np
import SimpleITK as sitk
import csv


def vectorise(path):
    arr = sitk.GetArrayFromImage(sitk.ReadImage(path))
    return arr.flatten()


config_csv = sys.argv[1]
outfile = sys.argv[2]
mask_path = sys.argv[3]

flat_mask = vectorise(mask_path)

img_arrays = []

with open(config_csv, 'rb') as csvfile,  open(outfile, 'wb') as csvout:
    reader = csv.reader(csvfile, delimiter=',')
    header = [x.strip() for x in reader.next()]
    writer = csv.writer(csvout, delimiter=',')
    writer.writerow(header[: -1] + ['intensity'])

    for row in reader:
        row = [r.strip() for r in row]
        if len(row) < 1:
            break
        d1_array = vectorise(row[-1].strip())
        for intensity_value, mask_value in zip(d1_array, flat_mask):
            if mask_value == 0.0:
                intensity_value = 'NA'
            writer.writerow([x.strip() for x in row[: -1]] + [intensity_value])