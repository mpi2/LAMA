#!/usr/bin/env python

"""
Lifted from http://statsmodels.sourceforge.net/ipdirective/_modules/scikits/statsmodels/sandbox/stats/multicomp.html
"""
import sys
import numpy as np




def fdrcorrection(pvals, alpha=0.05):

    pvals = np.asarray(pvals)

    num_nans = np.isnan(pvals).size

    pvals_sortind = np.argsort(pvals)
    pvals_sorted = pvals[pvals_sortind]
    sortrevind = pvals_sortind.argsort()

    ecdffactor = ecdf(pvals_sorted)

    pvals_corrected_raw = pvals_sorted / ecdffactor

    pvals_corrected_raw[np.isnan(pvals_corrected_raw)] = 100

    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]


    pvals_corrected[pvals_corrected>1] = 1
    return pvals_corrected[sortrevind]




def ecdf(x, numnans=0):
    '''no frills empirical cdf used in fdrcorrection
    '''
    nobs = len(x) + numnans
    return np.arange(1,nobs+1)/float(nobs)


infile = sys.argv[1]
outfile = sys.argv[2]

pvals = np.genfromtxt(infile, dtype=np.float128)


for q in fdrcorrection(pvals):
    print q