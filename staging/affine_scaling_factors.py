#!/usr/bin/env python


import utilities.transformations as trans
import numpy as np
import os
from os.path import join


TFORM_FILE_NAME = 'TransformParameters.0.txt'
AFFINE_INDENTIFIER = '(Transform "AffineTransform")'
SIMILARITY_INDENTIFIER = '(Transform "SimilarityTransform")'
TFORM_PARM_LINE_START = '(TransformParameters'


def get_scaling_factor(tform_params):

    if len(tform_parms) == 7: # Similarity
       return tform_parms[-1]

    mat = np.array(tform_params).reshape(4, 3).T
    hom = np.array([0, 0, 0, 1])
    hom_mat = np.vstack((mat, hom))
    s = trans.decompose_matrix(hom_mat)
    scaling = s[0][:3]

    return scaling[2]

def extract_affine_transformation_parameters(path):
    for file_ in os.listdir(path):
        if not file_ == TFORM_FILE_NAME:
            continue
        with open(os.path.join(path, file_)) as reader:
            first_line = reader.next()
            if not first_line.startswith((AFFINE_INDENTIFIER, SIMILARITY_INDENTIFIER)):
                continue
            for line in reader:
                if line.startswith(TFORM_PARM_LINE_START):
                    tform_params = line.strip().strip(')').split()[1:]
                    return [float(x) for x in tform_params]


if __name__ == '__main__':
    import sys
    folder = sys.argv[1]
    for sub_folder in os.listdir(folder):
        if os.path.isdir(join(folder, sub_folder)):
            tform_parms = extract_affine_transformation_parameters(os.path.join(folder, sub_folder))
            if tform_parms:
                scale_factor = get_scaling_factor(tform_parms)
                print "{}, {}".format(sub_folder, scale_factor)