#!/usr/bin/env python


import numpy as np
import SimpleITK as sitk
import pycircstat as pc
import os
from os.path import join, basename, splitext


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

wt = ([-29.54894447,  -5.97417927,  -0.65787357],
[-5.62563133,  0.25041604, -4.68448496],
[-6.81777334, -9.70918369, -8.01682091],
[-47.07903671,   5.7998538,   12.89711857],
[  1.58070254, -28.59241295, -18.62919044],
[-21.91890526,  10.54876423,  10.04658604],
[-11.60393143,  -9.75127697,  -7.42155218],
[ -0.55180448, -18.10203934, -22.81825638],
[ 6.25147772, -10.75350094,  -4.41338158])



mut = ( [ -41.53026581,  -9.98734188, -15.41089916],
        [ -41.53026581,  -9.98734188, -15.41089916],
        [ -41.53026581,  -9.98734188, -15.41089916])



# calculate the angel between
mut_angles = []
wt_angles = []

for wa in wt:
    angle = np.rad2deg(np.arctan2(wa[0], wa[1]))
    wt_angles.append(str(angle))

for ma in mut:
    angle = np.rad2deg(np.arctan2(ma[0], ma[1]))
    mut_angles.append(str(angle))

print "wt -> c({})".format(','.join(wt_angles))
print "mut -> c({})".format(','.join(mut_angles))
