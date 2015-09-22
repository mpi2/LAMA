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
import yaml


def vectorise(path):
    arr = sitk.GetArrayFromImage(sitk.ReadImage(path))
    return arr.flatten()

config_csv = sys.argv[1]
outfile = sys.argv[2]

with open(config_csv, 'rb') as csvfile,  open(outfile, 'wb') as csvout:
    reader = csv.reader(csvfile, delimiter=',')
    header = [x.strip() for x in reader.next()]
    writer = csv.writer(csvout, delimiter=',')
    writer.writerow(header[: -1] + ['intensity'])

    for row in reader:
        d1_array = vectorise(row[-1].strip())
        for value in d1_array:
            writer.writerow([x.strip() for x in row[: -1]] + [value])

