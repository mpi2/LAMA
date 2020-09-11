"""
usage

pytest -v -m "not notest" test_staging.py
"""
from pathlib import Path
from lama.staging import staging_metric_maker
from . import wt_registration_dir
import pytest


# @pytest.mark.notest
def test_otsu():
    """
    Test the creation of staging metrics using otsu thresholding to generate whole embryo masks

    Note: test data shoul dbe gernated frm lama in lama/tests/test_data/mutant_and_baseline_data/baseline
    This can be gerated by running the tests in lama/tests/test_data_generation.py

    """
    #original_vol_dir: Path, outdir: Path, otsudir: Path

    # Chose a WT reg dir to test
    rigid_reg_dir = wt_registration_dir / 'output/baseline/20140430_KLF7_E14.5_16.5f_WT_XY_rec_scaled_2.8248_pixel_14/output/registrations/rigid'
    outdir = rigid_reg_dir.parent.parent
    otsu_dir = outdir / 'otsu'
    staging_metric_maker.otsu_staging(rigid_reg_dir, outdir, otsu_dir)