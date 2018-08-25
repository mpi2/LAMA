#!/usr/bin/python

from . import hill_climber_c
import sys

label_map = sys.argv[1]
label_num = int(sys.argv[2])
jac_value = float(sys.argv[3])
padding = int(sys.argv[4])
out_dir = sys.argv[5]
ngen = sys.argv[6]


hill_climber_c.HcShrink(label_map, label_num, jac_value, padding, out_dir, ngen)