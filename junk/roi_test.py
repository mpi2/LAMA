#!/usr/bin/env python

from os.path import join
import os
import sys
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
from roi_overlay import make_normalization_roi_qc_images


roi = [[610, 204, 233], [620, 214, 243]]

outdir = '/home/neil/sig/LAMA_results/E14.5/120716_E14.5_14um_test_set/output/qc/normalisation_roi_overlays'
indir = '/home/neil/sig/LAMA_results/E14.5/120716_E14.5_14um_test_set/output/normalised_intensity_images'

make_normalization_roi_qc_images(indir, roi, outdir)