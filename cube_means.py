#!/usr/bin/env python

"""
Generate a statistics data file that can be used as input for vpv.py
"""
from simplejson import ordered_dict

import numpy as np
import SimpleITK as sitk
import argparse
import cPickle as pickle
import harwellimglib as hil
import sys
import os
from multiprocessing import Pool
import pprint
from collections import OrderedDict
import tempfile
from scipy import stats




def cube_vect_magnitude_mean(cube_of_vectors):
    """
    For a cube of deformation vectors, get the mean magnitude
    :param cube_of_vectors:
    :return: mean magnitude
    """
    vectors = []
    #Append each vector from the cube to a list
    for z in cube_of_vectors:
        for y in z:
            for vec in y:
                vectors.append(vec)
    #Get the mean vector. Then get the magnitude of it using np.linalg.norm
    return np.linalg.norm(np.mean(vectors, axis=0))


def cube_scalar_mean(cube_of_scalars):
    """
    For a cube of jacobian scalars, get the mean value
    :param cube_of_jac_scalars:
    :return:mean jacobian
    """
    return np.mean(cube_of_scalars)



def group_position_means(all_files_data):
    """
    Ggroup together mean magnitudes from each file by corresponding position.
    :param all_files_data: Data for all files at one cube position
    :return:pos_means, list
    """
    pos = [all_files_data[0][0]]  #=(x,y,z)
    all_means = []

    #add the each files cube average
    for deform_file_data in all_files_data:
        all_means.append(deform_file_data[1])

    pos.append(all_means)
    return pos


def run(WTs, mutants, analysis_type, chunksize, outfile):
    """
    :param jacobians:List of jacobian files
    :param deforms: List of defomation files
    :param outfile: Where to put the output
    :param threads: Max num of threads to use
    :param deform_file_pattern: defomation file regex finder in case these files are mixed in with others
    """

    if os.path.isdir(args.outfile):
        sys.exit("Supply an output file path not a directory")

    print('processing')

    wt_paths = hil.GetFilePaths(os.path.abspath(WTs))
    mut_paths = hil.GetFilePaths(os.path.abspath(mutants))

    if len(wt_paths) < 1:
        raise IOError("can't find volumes in {}".format(WTs))
    if len(mut_paths) < 1:
        raise IOError("can't find volumes in {}".format(mutants))

    print('### Wild types to process ###')
    print([os.path.basename(x) for x in wt_paths])
    print('### Mutants to process ###')
    print([os.path.basename(x) for x in mut_paths])

    vol_stats(wt_paths, mut_paths, analysis_type, chunksize, outfile)



def vol_stats(wts, muts, analysis_type, chunksize, outfile):

    memmapped_wts = memory_map_volumes(wts)
    memmapped_muts = memory_map_volumes(muts)

    # filename = os.path.basename(rawdata_file)
    # img = sitk.ReadImage(rawdata_file)
    # if data_type == 'intensities':
    #     print('normalizing')
    #     img = sitk.Normalize(img)
    zdim, ydim, xdim = memmapped_wts[0].shape[0:3]  # For vectors there's and extra dimension so can't just unpack

    # Create an array to store the
    ttest_result_array = np.zeros(shape=memmapped_wts[0].shape, dtype=np.float32)

    path_ = '/home/neil/share/registration_projects/120315_NXN/NXN_mutants/out/deformable/NXN_K1029-1_KO.mnc/NXN_K1029-1_KO.mnc.nii'

    test_array = sitk.GetArrayFromImage(sitk.ReadImage(path_))


    for z in range(0, zdim - chunksize, chunksize):
        for y in range(0, ydim - chunksize, chunksize):
            for x in range(0, xdim - chunksize, chunksize):


                wt_means = get_mean_cube(memmapped_wts, z, y, x, chunksize, analysis_type)
                mut_means = get_mean_cube(memmapped_muts, z, y, x, chunksize, analysis_type)

                ttest_result = ttest(wt_means, mut_means)

                #print ttest_result
                ttest_result_array[z:z+chunksize, y:y+chunksize, x:x+chunksize] = ttest_result

    result_vol = sitk.GetImageFromArray(ttest_result_array)
    sitk.WriteImage(result_vol, outfile)

def get_mean_cube(arrays, z, y, x, chunksize, a_type):

    means = []
    for arr in arrays:

        if a_type == 'def':
            cube_mean = cube_vect_magnitude_mean(arr[z:z+chunksize, y:y+chunksize, x:x+chunksize])
        else:
            cube_mean = cube_scalar_mean(arr[z:z+chunksize, y:y+chunksize, x:x+chunksize])
        means.append(cube_mean)

    return means


def ttest(wt, mut):
    """
    :param wt:
    :param mut:
    :return: float, pvalue
    """
    #Can't get scipy working on Idaho at the moment so use ttest from cogent package for now
    #return np.mean(wt) - np.mean(mut)
    return stats.ttest_ind(wt, mut)[0]

    #Trying out the Kolmogorov-Smirnov statistic on 2 samples.
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ks_2samp.html
    #return stats.ks_2samp(mut,  wt)[1]  # [1] is the pvalue

def memory_map_volumes(vol_paths):
    """
    Create memory-mapped volumes
    """
    memmapped_vols = []
    for vp in vol_paths:
        img = sitk.ReadImage(vp)
        # if norm:
        #     img = sitk.Normalize(img)
        array = sitk.GetArrayFromImage(img)
        tempraw = tempfile.TemporaryFile(mode='wb+')
        array.tofile(tempraw)
        memmpraw = np.memmap(tempraw, dtype=np.float32, mode='r', shape=array.shape)

        memmapped_vols.append(memmpraw)

    return memmapped_vols



if __name__ == '__main__':

    # mpl = multiprocessing.log_to_stderr()
    # mpl.setLevel(multiprocessing.SUBDEBUG)

    parser = argparse.ArgumentParser("messageS")
    parser.add_argument('-c', dest='cubesize', type=int, help='Size of the sub array', required=True)
    parser.add_argument('-w', dest='wt_vols_dir', help='Folder containing WT data ', required=True)
    parser.add_argument('-m', dest='mut_vols_dir', help='Folder containing Mut data', required=True)
    parser.add_argument('-a', dest='analysis_type', help='<int, def, jac> Intensities, deformation fields, or spatial jacobians', required=True)
    parser.add_argument('-o', dest='outfile', help='', required=True)
    # parser.add_argument('-r', dest='registered_vols', help='Folder containing registered vols, for intensity difference'
    #                                                        ' analysis', default=None)
    # parser.add_argument('-o', dest='outfile', help='File to store pickle file of means/stdvs of vectors,intensities etc', required=True)
    # parser.add_argument('-t', dest='threads', type=int, help='How many threads to use', default=4)
    # parser.add_argument('-dp', dest='deform_file_pattern', help='String that is contained in the deform file names',
    #                     default='deformationField')
    # parser.add_argument('-jp', dest='jac_file_pattern', help='String that is contained in the jacobian file names',
    #                     default='spatialJacobian')

    args = parser.parse_args()

    run(args.wt_vols_dir, args.mut_vols_dir, args.analysis_type, args.cubesize, args.outfile)
