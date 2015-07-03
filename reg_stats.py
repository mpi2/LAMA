#!/usr/bin/env python

"""
Generate a statistics data file that can be used as input for vpv.py

TODO: remove mask boilerplate
"""

import sys
import os
import argparse
import tempfile
import logging

import numpy as np
import SimpleITK as sitk
from scipy import stats
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import rpy2.robjects as robj
rstats = importr('stats')
import yaml

import harwellimglib as hil
import common

# can't install h5py on idaho at the moment
try:
    import h5py
except ImportError:
    print 'warning: cannot import h5py. Minc files cannot be analysed'

LOG_FILE = '_stats.log'
LOG_MODE = logging.DEBUG
TSCORE_OUT_SUFFIX = '_tscore.nrrd'


def reg_stats(config_path):
    """

    """

    # logfile = os.path.join(outfile, LOG_FILE)
    # common.init_log(os.path.dirname(logfile, 'Stats log', LOG_MODE))
    # common.init_log('Stats processing started')

    print('processing')

    # All paths in config file are relative to this
    config_dir = os.path.dirname(config_path)
    config = yaml.load(open(config_path, 'r'))
    n1 = config.get('n1')  # Do one against many analysis?
    mask = os.path.join(config_dir, config.get('fixed_mask'))

    stats_outdir = os.path.join(config_dir, 'stats')
    common.mkdir_force(stats_outdir)

    for analysis_name, reg_data in config['data'].iteritems():
        _analyse(analysis_name, reg_data, config_dir, stats_outdir, mask, n1)


def _analyse(analysis_name, reg_data, config_dir, stats_outdir, mask, n1):

    wt_dir = reg_data['wt']
    mut_dir = reg_data['mut']
    data_type = reg_data['datatype']  # Scalar/vector

    try:
        wt_img_paths = hil.GetFilePaths(os.path.join(config_dir, wt_dir))
        mut_img_paths = hil.GetFilePaths(os.path.join(config_dir, mut_dir))
    except OSError:
        sys.exit('Cant find or access the volumes')  # a bit extreme?

    if len(wt_img_paths) < 1:
        raise IOError("can't find volumes in {}".format(wt_dir))
    if len(mut_img_paths) < 1:
        raise IOError("can't find volumes in {}".format(mut_dir))

    print('### Wild types to process ###')
    print([os.path.basename(x) for x in wt_img_paths])
    print('### Mutants to process ###')
    print([os.path.basename(x) for x in mut_img_paths])

    analysis_out_dir = os.path.join(stats_outdir, analysis_name)
    common.mkdir_force(analysis_out_dir)

    many_against_many(wt_img_paths, mut_img_paths, data_type, analysis_out_dir, mask)
    if n1:
        one_against_many(wt_img_paths, mut_img_paths, data_type, analysis_out_dir, mask)


def one_against_many(wts, muts, data_type, analysis_dir, mask, memmap=False):

    blurred_wts = memory_map_volumes(wts, memmap, data_type)

    wt_stdvs = get_std(blurred_wts)

    for mut_path in muts:
        blurred_mut = memory_map_volumes([mut_path], memmap, data_type)[0]

        # Filter out any values below 2 standard Deviations
        blurred_mut[blurred_mut < wt_stdvs * 2] = 0

        mask_array(blurred_mut, mask, replace=0)
        img = sitk.GetImageFromArray(blurred_mut)
        out = os.path.join(analysis_dir, os.path.basename(mut_path))
        sitk.WriteImage(img, out)


def many_against_many(wts, muts, data_type, analysis_dir, mask, memmap=False):

    blurred_wts = memory_map_volumes(wts, memmap, data_type)
    blurred_muts = memory_map_volumes(muts, memmap, data_type)

    try:
        mask_img = sitk.ReadImage(mask)
        mask_arr = sitk.GetArrayFromImage(mask_img)
    except RuntimeError:
        mask_arr = None

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
    outfile = os.path.join(analysis_dir, TSCORE_OUT_SUFFIX)

    sitk.WriteImage(t_img, outfile)


def mask_array(array, mask, replace=0):
    """
    mask a ndarray

    parameters
    ----------
    array: ndarray
        array to be masked
    mask: SimpleITK Image
        Image mask
    replace: Anything that can fit in a numpy array

    Returns
    -------
    Modifies ndarray in place
    """
    mask_img = sitk.ReadImage(mask)
    mask_arr = sitk.GetArrayFromImage(mask_img)
    array[mask_arr == 0] = replace


def get_std(arrays):
    """
    Return an ndarray of standard deviations derived element-wise from a series of input arrays

    Parameters
    ----------
    arrays: array
        array of ndarrays

    Returns
    -------
    ndarray
        array of standar deviations
    """
    stdv = np.std(arrays, axis=0)
    return stdv


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
    """

    """
    if type(mask) == np.ndarray:
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
    if analysis_type == 'vector':
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

    parser = argparse.ArgumentParser("Stats component of the phenotype detection pipeline")
    parser.add_argument('-c', '--config', dest='config', help='yaml config file', required=True)

    args = parser.parse_args()

    reg_stats(args.config)
