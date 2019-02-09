#!/usr/bin/env python

"""
Henrik: Reimplementing this to be a hill climber
"""



import SimpleITK as sitk

def test_DJD(label_path, ngen):

    im = sitk.ReadImage(label_path)

    for i in range(ngen):
        jacobian = sitk.DisplacementFieldJacobianDeterminant(im)


if __name__ == '__main__':
    import sys

    label = sys.argv[1]
    ngen = int(sys.argv[2])

    test_DJD(label, ngen)



