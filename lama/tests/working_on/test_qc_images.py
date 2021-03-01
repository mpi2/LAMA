"""
Functions to test QC image creation.
Manual checking of the output is required to make sure the images are as expected.


Notes
-----
 These tests require test_data_generation.py to have been run previously so it has a data folder to work on
 The test data is very small, but it should give an idea of whether the QC images are being made correclty
 I may add full size images in the future.
"""

import pytest
from lama.tests import wt_registration_dir, target_img
from lama.qc import qc_images

# A lama specimen directory created when running test_data_generation.py
lama_specimen_dir = wt_registration_dir / 'output' / 'baseline' / '20140122_SCN4A_18.1_e_wt_rec_scaled_3.1241_pixel_14'
outdir = lama_specimen_dir / 'output' / 'qc'

def test_qc_images():
    qc_images.make_qc_images(lama_specimen_dir, target_img, outdir)
    print('p')