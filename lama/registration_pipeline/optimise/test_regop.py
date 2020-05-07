
from pathlib import Path
import lama
from . import regop_points
import pytest
import shutil
import os
import pytest


test_dir = Path(lama.__file__).parent / 'tests' / 'test_data' / 'points_test_data'

# Note a fully test directory needs to be downloaded to here: TODO: put that online somewhere
reg_dir = test_dir / '20131025_CBX2_E14.5_12.1g_WT_XY_rec_scaled_6.8627_pixel_13'

target_points = reg_dir / 'popavg_points.fcsv'
elx_points = reg_dir / 'popavg_points.txt'
trformed_elx_ponts =  reg_dir / 'popavg_points_transformed.txt'
moving_points = reg_dir / 'moving_points.fcsv'


@pytest.mark.skip
def test_convert_points():

    shutil.rmtree(elx_points, ignore_errors=True)

    regop_points.slicer_points_to_elastix(target_points, elx_points)

    assert elx_points.is_file()


# @pytest.mark.skip
def test_transform_point():

    new_points_file = regop_points._transform_points(reg_dir, elx_points)
    assert new_points_file.is_file()

@pytest.mark.skip
def test_compare_points():
    target_inverted_points = reg_dir / 'output' / 'inverted_points' / 'outputpoints.txt'
    regop_points._compare_points(target_inverted_points, moving_points)

