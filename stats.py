#!/usr/bin/env python

"""
Generate a statistics data file that can be used as input for vpv.py
"""


import numpy as np
import SimpleITK as sitk
import argparse
import harwellimglib as hil
import sys
import os
import tempfile
from scipy import stats
import math
import h5py
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
rstats = importr('stats')





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
    shape = memmapped_wts[0].shape[0:3]  # For vectors there's and extra dimension so can't just unpack
    dtype = memmapped_wts[0].dtype
    zdim, ydim, xdim = shape

    print("Calculating statistics")

    #pvalues = np.zeros(shape, dtype=memmapped_wts[0].dtype)
    tstats = []
    pvalues = []

    for z in range(0, zdim - chunksize, chunksize):
        for y in range(0, ydim - chunksize, chunksize):
            for x in range(0, xdim - chunksize, chunksize):

                wt_means = get_mean_cube(memmapped_wts, z, y, x, chunksize, analysis_type)
                mut_means = get_mean_cube(memmapped_muts, z, y, x, chunksize, analysis_type)
                tstat, pval = ttest(wt_means, mut_means)

                pvalues.append(pval)
                tstats.append(tstat)

                # if not np.isfinite(tstat):
                #     print tstat, pval, wt_means, mut_means
                #     sys.exit()

    print("calculating FDR")
    qvalues = fdr(pvalues)

    # Filter out t statistics that have corresponding pvalues > 0.
    filtered_tstats = filter_tsats(tstats, qvalues)

    no_inf = remove_infinite(filtered_tstats)

    tstat_volume = expand(no_inf, shape, chunksize)

    result_tstats = sitk.GetImageFromArray(tstat_volume)
    sitk.WriteImage(result_tstats, outfile)


def remove_infinite(tstats):
    t = np.array(tstats)
    pos_mask = (~np.isfinite(t)) & (t > 0)
    neg_mask = (~np.isfinite(t)) & (t < 0)
    t[pos_mask] = 0
    t[neg_mask] = 0
    t[pos_mask] = t.max()
    t[neg_mask] = t.min()
    return t

def fdr(pvalues):

    qvals = rstats.p_adjust(FloatVector(pvalues), method='BH')
    return qvals


def expand(qstats, shape, chunksize):
    tstats_count = 0
    expanded = np.zeros(shape, dtype=np.float)

    for z in range(0, shape[0] - chunksize, chunksize):
        for y in range(0, shape[1] - chunksize, chunksize):
            for x in range(0, shape[2] - chunksize, chunksize):
                expanded[z:z+chunksize, y:y+chunksize, x:x+chunksize] = qstats[tstats_count]
                tstats_count += 1

    return expanded


def filter_tsats(tstats, qvalues):
    """
    Convert to numpy arrays and set to zero any tscore that has a corresponding pvalue > 0.05
    :param tsats: array
    :param qvalues: array
    :return: np.ndarray, filtered t statistics
    """
    mask = qvalues > 0.05
    tstats[mask] = 0
    return tstats



def get_mean_cube(arrays, z, y, x, chunksize, a_type):
    """

    :param arrays:
    :param z:
    :param y:
    :param x:
    :param chunksize:
    :param a_type:
    :return:
    """
    means = []
    for arr in arrays:

        if a_type == 'def':
            means.append(np.mean(arr[z:z+chunksize, y:y+chunksize, x:x+chunksize]))
        else:
            means.append(np.mean(arr[z:z+chunksize, y:y+chunksize, x:x+chunksize]))
    return means



def ttest_bu(wt, mut):

    wilcox = rstats.wilcox_test(FloatVector(wt), FloatVector(mut))
    try:
        statistic = wilcox[0]
        p_value = wilcox[2]
    except IndexError:
        # If test fails (eg both vectors are all zeros
        statistic = 0.0
        p_value = 1
    return statistic, p_value
    # stat, pval = stats.ttest_ind(wt, mut)
    # if math.isnan(pval):
    #     pval = 1.0
    # return stat, pval





def ttest(wt, mut):
    """
    Calculate the pvalue and the tstatistic for the wt and mut subsampled region

    Args:
       wt (list):  wildtype values
       mut (list): mutant values

    Returns:
       tuple: (pvalue(float), is_tscore_positive(bool)
    """
    #tscore, pval = stats.ttest_ind(wt, mut)[0:2]
    tscore, pval = stats.ttest_ind(wt, mut)[0:2]

    # set the pvalue to negative if the t-statistic is negative
    # if tscore < 0:
    #     is_tscore_positive = False
    # else:
    #     is_tscore_positive = True

    # Set pval nan values to 1. This can happen if all values in tested arrays are 0
    if math.isnan(pval):
        pval = 1.0

    return tscore, pval


def memory_map_volumes(vol_paths):
    """
    Create memory-mapped volumes
    """
    memmapped_vols = []
    for vp in vol_paths:
        if vp.lower().endswith('mnc'):
            img = read_mnc(vp)
        else:
            img = sitk.ReadImage(vp)
        np_array = sitk.GetArrayFromImage(img)
        data_type = np_array.dtype

        array = sitk.GetArrayFromImage(img)
        tempraw = tempfile.TemporaryFile(mode='wb+')
        array.tofile(tempraw)
        memmpraw = np.memmap(tempraw, data_type, mode='r', shape=array.shape)

        memmapped_vols.append(memmpraw)

    return memmapped_vols

def read_mnc(minc_path):
    """
    Read in a mnc file
    Args:
        minc_path (str): path to mnc file
    Return:
        numpy array
    """
    minc = h5py.File(minc_path, "r")['minc-2.0']
    nparray = minc['image']['0']['image']
    sitk_array = sitk.GetImageFromArray(nparray)
    return sitk_array





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
