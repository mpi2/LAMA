#!/usr/bin/python


from utilities import glcm3d
import reg_stats
import numpy as np

def run(wt, mut, mask, out):
    reg_stats.calculate_glcm_metrics(wt, mut, mask, out)



if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser("testing glcm paramters")
    parser.add_argument('-w', dest='wt', help='Numpy file of wild type glcms', required=True)
    parser.add_argument('-m', dest='mut', help='Numpy file of mutant type glcms', required=True)
    parser.add_argument('-mask', dest='mask', help='Mask volume path', required=True)
    parser.add_argument('-o', dest='out', help='out dir', required=True)

    args = parser.parse_args()

    run(args.wt, args.mut, args.mask, args.out)
