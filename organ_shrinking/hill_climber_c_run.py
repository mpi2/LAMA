#!/usr/bin/python

import hill_climber_c
import sys

label = sys.argv[1]
jac_value = sys.argv[2]
out_dir = sys.argv[3]
ngen = sys.argv[4]


hill_climber_c.HcShrink(label, jac_value, out_dir, ngen)