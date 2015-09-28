"""
Given a subarray cube, and the orginal dimensions of that subarray, insert it back into a vactor image of the same
dimensiuons of the one it cam from
"""

import numpy as np
import SimpleITK as sitk

sub_array_dims = ((170, 205), (65, 91), (69, 99))
original_dims = [360, 256, 171]
sub_array_path = '/home/neil/work/defined_abnormalites/270915_single_thymus/shrink_0.7_270915/out/intermediate_results/def_0.nrrd'
out = '/home/neil/work/defined_abnormalites/270915_single_thymus/shrink_0.7_270915/out/intermediate_results/deformation.nrrd'


sub_array = sitk.GetArrayFromImage(sitk.ReadImage(sub_array_path))

def_field = np.zeros(original_dims + [3])





