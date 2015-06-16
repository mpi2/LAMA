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
import rpy2.robjects as robj
rstats = importr('stats')
import subprocess


def reg_stats(wildtypes, mutants, analysis_type, outfile, mask, memmap=False):
    """
    :param jacobians:List of jacobian files
    :param deforms: List of defomation files
    :param outfile: Where to put the output
    :param threads: Max num of threads to use
    :param deform_file_pattern: defomation file regex finder in case these files are mixed in with others
    """
    log_path = os.path.join(os.path.split(outfile)[0], 'stats.log')
    log = open(log_path, 'w', 0)
    script_base = os.path.split(os.path.realpath(__file__))[0]
    git_path = os.path.join(script_base, 'git')

    try:
        git_version = subprocess.check_output(['git', "rev-parse", 'HEAD'], cwd=git_path)
    except (subprocess.CalledProcessError, OSError):
        git_version = "stats version not available"

    log.write("Harwell anotomical phenotype detection pipeline\n"
              "Stats module - version: {}".format(git_version))

    if os.path.isdir(outfile):
        sys.exit("Supply an output file path not a directory")

    print('processing')

    try:
        wt_paths = hil.GetFilePaths(os.path.abspath(wildtypes))
        mut_paths = hil.GetFilePaths(os.path.abspath(mutants))
    except OSError:
        sys.exit('Cant find or access the volumes')

    if len(wt_paths) < 1:
        raise IOError("can't find volumes in {}".format(wildtypes))
    if len(mut_paths) < 1:
        raise IOError("can't find volumes in {}".format(mutants))

    print('### Wild types to process ###')
    print([os.path.basename(x) for x in wt_paths])
    print('### Mutants to process ###')
    print([os.path.basename(x) for x in mut_paths])

    vol_stats(wt_paths, mut_paths, analysis_type, outfile, mask, memmap)



def vol_stats(wts, muts, analysis_type, outfile, mask, memmap=False):

    blurred_wts = memory_map_volumes(wts, memmap, analysis_type)
    blurred_muts = memory_map_volumes(muts, memmap, analysis_type)
    if mask:
        mask_img = sitk.ReadImage(mask)
        mask_arr = sitk.GetArrayFromImage(mask_img)

    shape = blurred_wts[0].shape[0:3]  # For vectors there's and extra dimension so can't just unpack

    print("Calculating statistics")
    tstats, pvalues = ttest(blurred_wts, blurred_muts)

    print("Calculating FDR")
    qvalues = fdr(pvalues, mask_arr)

    print 'q', min(qvalues), max(qvalues), np.mean(qvalues)
    # reshape

    filtered_t = filter_tsats(tstats, qvalues)

    t_vol = filtered_t.reshape(shape)

    t_img = sitk.GetImageFromArray(t_vol)

    sitk.WriteImage(t_img, outfile)


def get_vector_magnitudes(img):
    """
    For a cube of deformation vectors, get the mean magnitude
    :param cube_of_vectors:
    :return: mean magnitude
    """
    print "getting deformation magnitudes"
    arr = sitk.GetArrayFromImage(img)
    #Get the mean vector. Then get the magnitude of it using np.linalg.norm
    scalars = np.sqrt((arr*arr).sum(axis=3))
    return sitk.GetImageFromArray(scalars)


def fdr(pvalues, mask):

    flat_mask = mask.flatten()
    pvalues[flat_mask == 0] = robj.NA_Real
    qvals = np.array(rstats.p_adjust(FloatVector(pvalues), method='BH'))
    qvals[np.isnan(qvals)] = 1
    return qvals


def filter_tsats(tstats, qvalues):
    """
    Convert to numpy arrays and set to zero any tscore that has a corresponding pvalue > 0.05
    :param tsats: array
    :param qvalues: array
    :return: np.ndarray, filtered t statistics
    """
    assert len(tstats) == len(qvalues)
    t = np.array(tstats)
    q = np.array(qvalues)
    mask = q > 0.07
    t[mask] = 0

    return t


def ttest(wt, mut):
    """
    Calculate the pvalue and the tstatistic for the wt and mut subsampled region

    Args:
       wt (list):  wildtype values
       mut (list): mutant values

    Returns:
       tuple: (pvalue(float), is_tscore_positive(bool)
    """
    wt_flat = flatten(wt)
    mut_flat = flatten(mut)

    tscores, pvals = stats.ttest_ind(mut_flat, wt_flat)

    mask_t = np.isnan(tscores)
    tscores[mask_t] = 0
    mask_p = np.isnan(pvals)
    pvals[mask_p] = 1

    return tscores, pvals

def flatten(arrays):
    one_d = []
    for arr in arrays:
        f = arr.flatten()
        one_d.append(f)
    stacked = np.vstack(one_d)
    return stacked


def memory_map_volumes(vol_paths, memmap=False, analysis_type='int'):
    """
    Create memory-mapped volumes
    """

    vols = []
    for vp in vol_paths:
        img = sitk.ReadImage(vp)
        blurred = blur(img, analysis_type)
        array = sitk.GetArrayFromImage(blurred)
        data_type = array.dtype
        if memmap:
            tempraw = tempfile.TemporaryFile(mode='wb+')
            array.tofile(tempraw)
            memmpraw = np.memmap(tempraw, data_type, mode='r', shape=array.shape)

            vols.append(memmpraw)
        else:
            vols.append(array)
    return vols


def blur(img, analysis_type):
    if analysis_type == 'def':
        img = get_vector_magnitudes(img)
    # previous: 1.0, 8, 0.001
    blurred = sitk.DiscreteGaussian(img, 0.5, 4, 0.01, False)
    return blurred

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
    parser.add_argument('-w', dest='wt_vols_dir', help='Folder containing WT data ', required=True)
    parser.add_argument('-m', dest='mut_vols_dir', help='Folder containing Mut data', required=True)
    parser.add_argument('-mask', dest='mask', help='Population average mask', required=False)
    parser.add_argument('-a', dest='analysis_type', help='<int, def, jac> Intensities, deformation fields, or spatial jacobians', required=True)
    parser.add_argument('-o', dest='outfile', help='', required=True)
    parser.add_argument('-mem', dest='memmap', help='', action="store_true", default=False)

    args = parser.parse_args()



    reg_stats(args.wt_vols_dir, args.mut_vols_dir, args.analysis_type, args.outfile, args.mask, args.memmap)
