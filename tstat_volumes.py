#!/usr/bin/env python
"""
This module reads in a previosly computed file that contains the subsampled
mean displacemnt vector, jacobian, or intensity data from the wild type registration. Each

Will output
"""

import pprint
import argparse
import pickle
import sys
from scipy import stats
import os
import numpy as np
import SimpleITK as sitk
import harwellimglib as hil
from collections import namedtuple
import pprint

def volume_compare_dirs(wt_dir, mut_dir, out_dir, )


def volume_compare(mut_deform_stats, wt_deform_stats, outdir, mut_deform_dir, filename_pattern, pvalue, data_type):
    """
    :param mut_deform_stats:
    :param wt_deform_stats:
    :param deform_vectors_dir:
    :param mut_deforms:
    :param wt_data:
    :param outdir:
    :return:
    """

    chunksize, wt_data, mut_data = _get_data(mut_deform_stats, wt_deform_stats)
    pos_pvals, deform_dimensions = _compare(chunksize, wt_data, mut_data)

    mean_deform = 'notused'  # Need to get deformation fields only if looking at vectors. Can't do truth on ndarray
    if data_type == "deformation":
        deform_paths = hil.GetFilePaths(mut_deform_dir, pattern=filename_pattern)
        if len(deform_paths) < 1:
            sys.exit("No deformation files found in folder containg {} in filename".format(filename_pattern))
        mean_deform = _mean_deform_vols(deform_paths)

    _output_volume(pos_pvals, outdir, chunksize, deform_dimensions, pvalue, mean_deform, data_type)


def _get_data(mut_deform_stats, wt_deform_stats):
    """
    Unpickle the the cube data (need to find a new name for thesse data) and do some checking
    :param mutant_mean_cubes: dict {mutid:[cubes]....}
    :param savedwt_data: file path
    :return:
    """

    with open(wt_deform_stats, 'rb') as wt_pickle, open(mut_deform_stats, 'rb') as mut_pickle:
        wt_data = pickle.load(wt_pickle)
        mut_data = pickle.load(mut_pickle)

    if wt_data['data_type'] != mut_data['data_type']:
        sys.exit("Looks like you are comparing different data types. Jacobians and deformation vectors maybe")

    if wt_data['chunksize'] != mut_data['chunksize']:
        sys.exit("the input data files do not have the same chunck size. Remake them")

    if wt_data['deform_dimensions'] != mut_data['deform_dimensions']:
        sys.exit("the input data files do not have the same dimensions. Remake them")

    chunksize = wt_data['chunksize']
    return chunksize, wt_data, mut_data


def _mean_deform_vols(paths):
    """
    Create a mean vector field file from a bunch of vector files
    :param vol_dir:
    :return:
    """
    first_path = paths.pop(0)
    vol = sitk.ReadImage(first_path)
    summed_array = sitk.GetArrayFromImage(vol)
    count = 1

    for path in paths:
        current_vol = sitk.ReadImage(path)
        current_array = sitk.GetArrayFromImage(current_vol)
        summed_array = summed_array + current_array
        count += 1

    mean_deform_array = summed_array / count

    return mean_deform_array


def _compare(chunksize, wt_data, mut_data):
    """
    loop over both sets of data and get some statistics for each cube
    :param chunksize:
    :param wt_data:
    :param mut_data:
    :return:
    """
    cube_results = []
    cubestats = namedtuple('Cubestats', ['pos', 'pval', 'mutmean', 'wtmean'])
    wt_means = wt_data['data']
    mut_means = mut_data['data']


    for (wt_index, wt_cube), (mut_index, mut_cube) in zip(wt_means.items(), mut_means.items()):
        #DEBUG: wt_cube and mut_cube botth == (0,0,0)
        assert wt_index == mut_index  # Make sure we are comapring the same positions
        pvalue = ttest(wt_cube, mut_cube)
        #print mut_cube, np.mean(mut_cube)
        cs = cubestats(wt_index, pval=pvalue, wtmean=np.mean(wt_cube), mutmean=np.mean(mut_cube))
        cube_results.append(cs)

    return cube_results, mut_data['deform_dimensions']


def _rgb(minimum, maximum, value):
    minimum, maximum = float(minimum), float(maximum)
    halfmax = (minimum + maximum) / 2
    r = int(max(0, 255*(1 - value/halfmax)))
    b = 255 - r
    g = 0
    return r, g, b


def _output_volume(pos_pvals_mean, outpath, chunksize, deform_dimensions, p_cuttoff, mean_deform, datatype):
    """

    :param pos_pvals_mean: Dict. {named_tuple}
    :param outpath: Where to put the results
    :param chunksize:
    :param deform_dimensions:
    :param p_cuttoff:
    :param mean_deform: Mean deformation field made from the input mutant. Is 'None' if looking at jacobians etc.
    :return:
    """
    #result array is a label volume. Remove a dimesion if working with vectors
    if len(deform_dimensions) > 3:
        deform_dimensions = deform_dimensions[0:3]
    result_array = np.zeros(deform_dimensions, dtype=np.uint8)

    passed_cube_list = []  # For outputing a list of readable passed cubes

    # The output pixel scale will go from (240,240,0) - yellow, to (240,0,0) - red
    # The approx range of pvals is 0.05 to 0.000005
    # old_min = 0.000005
    # old_max = p_cuttoff

    for cubestats in pos_pvals_mean:
        label = vox_val = 0
        pos = cubestats.pos
        pvalue = cubestats.pval

        if pvalue <= p_cuttoff:
            if datatype == 'jacobian':
                label = _jacobian_labelling(cubestats.mutmean)
            if datatype == 'intensities':
                label = intensity_labelling(cubestats.mutmean, cubestats.wtmean)
            if datatype == 'deformation':
                label = _deformation_labelling(cubestats.mutmean, cubestats.wtmean)

            passed_cube_list.append('pos: {0} pval: {1}, wtmean: {2}, mutmean: {3}'.format(
                pos, pvalue, cubestats.wtmean, cubestats.mutmean))

        #A label volume. 0: not significant, 1: contraction (jac< 1) 2: expansion (jac>1)
        result_array[
            pos[0]: pos[0] + chunksize,  # z extent
            pos[1]: pos[1] + chunksize,  # y extent
            pos[2]: pos[2] + chunksize   # x extent
        ] = label

        #------------------------------------------#
        #This bit is just for filtering out vectors
        #------------------------------------------#
        if mean_deform != 'notused':
            if not pvalue <= p_cuttoff:  # If the cube does not pass, set the vectors here to 0
                #A volume containing mean deform vectors that pass the pvalue cutoff
                mean_deform[
                    pos[0]: pos[0] + chunksize,  # z extent
                    pos[1]: pos[1] + chunksize,  # y extent
                    pos[2]: pos[2] + chunksize   # x extent
                ] = vox_val

    result_vol = sitk.GetImageFromArray(result_array)
    sitk.WriteImage(result_vol, os.path.join(outpath, 'seg_{0}_p{1}.tiff'.format(datatype, p_cuttoff)))

    if mean_deform != 'notused':
        filtered_deform = sitk.GetImageFromArray(mean_deform)
        sitk.WriteImage(filtered_deform, os.path.join(outpath, 'filtered_deformation_field.mhd'))

    passed_cube_list.insert(0, 'number of cubes pvalue < {0}: {1}'.format(p_cuttoff, len(passed_cube_list)))
    with open(os.path.join(outpath, 'hits_{0}_p{1}.txt'.format(datatype, p_cuttoff)), 'w') as fh:
        pprint.pprint(passed_cube_list, fh)


def _deformation_labelling(mutmean, wtmean):
    """
    Not sure about the best way to label yet. For now just stick with assigning labels by whether they are larger
    or smaller in the mutant
    :return:
    """
    if mutmean <= wtmean:
        label = 1
    if mutmean > wtmean:
        label = 2
    return label


def _jacobian_labelling(mut_mean_jac):
    """
    < 1 --> contraction
    > 1 --> expansion
    :param mut_mean_jac:
    :return:
    """
    if mut_mean_jac <= 1.0:
        label = 1  # Contratction
    else:
        label = 2  # Expansion
    return label


def intensity_labelling(mutmean, wtmean):
    """
    1 --> less intense in mutant
    2 --> more intense in mutant
    :param wtmean:
    :param mutmean:
    :return:
    """
    if mutmean < wtmean:
        label = 1
    if mutmean > wtmean:
        label = 2
    return label


def ttest(wt, mut):
    """
    :param wt:
    :param mut:
    :return: float, pvalue
    """
    #Can't get scipy working on Idaho at the moment so use ttest from cogent package for now
    return stats.ttest_ind(mut, wt)[1]

    #Trying out the Kolmogorov-Smirnov statistic on 2 samples.
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ks_2samp.html
    #return stats.ks_2samp(mut,  wt)[1]  # [1] is the pvalue

#TODO: Complain if proper datatype is not given. Better still have a flag specific for each datatype
if __name__ == '__main__':
    parser = argparse.ArgumentParser("Compare mutant displacement vectors to that of wild types")
    parser.add_argument('-m', dest='mut_deform_stats', help='mutant deformation file',
                        required=True)
    parser.add_argument('-w', dest='wt_deform_stats', help='location of the saved wild type data',
                        required=True)
    parser.add_argument('-o', dest='outdir', help='Dir to save outfile and volume',
                        required=True)
    parser.add_argument('-d', dest='deform_vector_dir', help='Dir with mutant deform vectors. Only needed for vector'
                                                             'analysis')
    parser.add_argument('-p', dest='pvalue', help='P value to filter vectors by', type=float,
                        default=0.05/5000)
    parser.add_argument('-fp', dest='filename_patternmatch', help='pattern to search for deform files',
                        default='deformationField')
    parser.add_argument('-dt', dest='data_type', help='What are you analysing? <vectors> <jacobians> <intensities>',
                        required=True)
    args = parser.parse_args()

    if args.data_type not in ['deformation', 'jacobian', 'intensities']:
        sys.exit("-dt needs to be one of <deformation> <jacobian> <intensities>")

    volume_compare(args.mut_deform_stats, args.wt_deform_stats, args.outdir, args.deform_vector_dir, args.filename_patternmatch,
        args.pvalue, args.data_type)




