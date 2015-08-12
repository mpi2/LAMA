#!/usr/bin/env python
from wx._core import wxEVT_MOVING

import SimpleITK as sitk
import numpy as np
import os
from nipype.interfaces.nipy.utils import Similarity
#from tempfile import NamedTemporaryFile

def run(fixed, moving, out):
    f = sitk.ReadImage(fixed)

    m = sitk.ReadImage(moving)
    moving = sitk.GetArrayFromImage(m)


    # rotate

    max_corr = {'turns': 0,
                'flip': False,
               'corr': 0}

    fixed_nii = "FixedTemp.nii"
    sitk.WriteImage(f, fixed_nii)

    tempMoving = "tempMoving.nii"

    for turn in [1, 2, 3, 4]:
        rotated = np.rot90(moving, turn)
        rotated_img = sitk.GetImageFromArray(rotated)
        sitk.WriteImage(rotated_img, tempMoving)

        corr = get_sim(fixed_nii, tempMoving)

        sitk.WriteImage(rotated_img, os.path.join(out, "{}.nrrd".format(str(turn))))
        print turn, corr
        if corr > max_corr['corr']:
            max_corr['corr'] = corr
            max_corr['turns'] = turn
            max_corr['flip'] = False

        flipped = np.flipud(rotated)

        flipped_img = sitk.GetImageFromArray(flipped)
        sitk.WriteImage(flipped_img, tempMoving)

        corr = get_sim(fixed_nii, tempMoving)

        sitk.WriteImage(flipped_img, os.path.join(out, "{}fliupped.nrrd".format(str(turn))))
        print turn, corr
        if corr > max_corr['corr']:
            max_corr['corr'] = corr
            max_corr['turns'] = turn
            max_corr['flip'] = True


    result_arr = np.rot90(moving, max_corr['turns'])
    if max_corr.get('flip'):
        result_arr = np.flipud(result_arr)

    result_img = sitk.GetImageFromArray(result_arr)
    sitk.WriteImage(result_img, os.path.join(out, "best_{}.nrrd".format(str(max_corr['turns']))))


def get_sim(fixed, moving):
    similarity = Similarity()
    similarity.inputs.volume1 = fixed
    similarity.inputs.volume2 = moving
    similarity.inputs.metric = "cc"
    res = similarity.run()



    return res.outputs.similarity



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("")
    parser.add_argument('-f', dest='fixed', help='fixed', required=True)
    parser.add_argument('-m', dest='moving', help='moving', required=True)
    parser.add_argument('-o', dest='out', help='out', required=True)
    args = parser.parse_args()

    run(args.fixed, args.moving, args.out)


