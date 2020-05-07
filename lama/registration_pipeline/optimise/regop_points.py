"""
The aim of this module is to compare distances from pairs of points to assess registration accuracy.

Two sets of points need to be defined:
    fixed points: points defined on your fixed volume (population average)
    moving points: points defined on your moving volume

Using a registrtion output of Lama, the transform parameters that defien the transform of moving -> fixed are
used to transform the fixed points back onto the input moving volume (on which the moving points were defind)

Note: Iversion of transform parameters is not required as elastix can transform fixed points to moving by using only
the forward transform parameter files. So to speed things up,  use "skip_transform_inversion: true" in the Lama config.

For each pair of points on the moving and transformed fixed points, the euclidean distance is computed, which are then
written to target_moving_distances.csv

The transformed fixed points are written to "transformed_target_points.fcsv" in Slicer fiducial format so they can be
loaded onto the moving image to get a visual clue to regisrtation accuracy.

If using as a module, run() is the only function that needs to be called. It also has a CLI.


Examples
--------
CLI
regop_points.py -d data/baselines/output/baseline/20190807_SCN1B-DEL567-EM1_E14.5_11.3c/ -f pop_average_points.fcsv -m 11_3c.fcsv

"""

from typing import List
from pathlib import Path
import subprocess as sub
import pandas as pd
from lama.registration_pipeline.optimise.points import Points
import os
import shutil
from scipy import spatial


def run(lama_reg_dir: Path, target_points: Path, moving_points: Path, out_dir: Path):
    """
    The main function in this module.

    Parameters
    ----------
    lama_reg_dir
        A Lama registration output directory. Minimal requirement is that it contains
         'output/registrations/', which must contain reg_order.txt that should created by lama

            example reg_order.txt
                rigid
                similarity
                affine
                deformable_192_to_10

    target_points
        A Slicer fiducial file of points defined on the target used in the registraiton in lama_reg_dir
    moving_points
        A Slicer fiducial file of points defined on the miving image that was registered in lama_reg_dir
    out_dir
        if not specified, output will go to lama_reg_dir/output/inverted_points and required dirs will be made
    """

    if not out_dir:
        out_dir = lama_reg_dir / 'output' / 'inverted_points'

    out_dir.mkdir(exist_ok=True)

    target_elastix_input = target_points.with_suffix('.txt')

    # Convert slicer target points to elastix input format
    Points().set_slicer(target_points).write_elastix_input(target_elastix_input)

    # Transform the target points
    tformed_elx_points = _transform_points(lama_reg_dir, target_elastix_input, out_dir)

    # Get the pairwise distances between points
    distances = _compare_points(tformed_elx_points, moving_points)
    distance_file = out_dir / 'target_moving_distances.csv'
    distances.to_csv(distance_file)

    print(f'distance results written too: {distance_file}')

    # Create Slicer file to look at the transformed points
    tformed_target_slicer_file = out_dir / 'transformed_target_points.fcsv'
    Points().set_elastix_output(tformed_elx_points).write_slicer(tformed_target_slicer_file)


def _transform_points(reg_dir: Path, points_file: Path, out_dir: Path) -> Path:
    """
    Given a Lama output regisrtation directory and a points file, apply the concatenated transforms to the points
    and return the output points file.

    Parameters
    ----------
    reg_dir
        Lama output directory
    points_file
        elastix points to be transformed
    out_dir
        if None create deafult folder with reg_dir

    Returns
    -------
    elastix transformed points
    """

    tforms = {}
    reg_order_file = reg_dir / 'output' / 'registrations' / 'reg_order.txt'
    reg_order = []

    with open(reg_order_file, 'r') as fh:
        for line in fh:
            if line.strip():
                stage = line.strip()
                reg_order.append(stage)
                tform = reg_dir / 'output' / 'registrations' / stage / reg_dir.name / 'TransformParameters.0.txt'
                tforms[stage] = tform

    first_tp_file = _modify_tforms(tforms, out_dir)

    os.chdir(out_dir)  # transformix will look for paths relative to cwd

    cmd = [
        'transformix',
        '-def', str(points_file),
        '-tp', str(first_tp_file),
        '-out', out_dir
    ]

    try:
        sub.check_output(cmd)
    except Exception as e:
        print('transformix failed with this command: {}\nerror message:'.format(cmd))
        raise

    return out_dir / 'outputpoints.txt'


def _modify_tforms(tforms: List[Path], outdir) -> Path:
    """
    Concatenate multiple elastix transforms by adding InitialTransformParameterFileName
    to point to previous TP files

    Returns
    -------
    The path of the final transform prameter file that links to all others.
    """

    # copy tp files to new folder
    new_names = {}
    for stage_id, tpfile in tforms.items():
        shutil.copyfile(tpfile, outdir / f'{stage_id}.txt')
        new_names[stage_id] = outdir / f'{stage_id}.txt'

    if len(tforms) < 2:  # Cannot have initial tform file as we need at least 2
        raise(ValueError('At least two transform parameter files needed'))

    for i, stage in enumerate(reversed(list(new_names.keys()))):
        tp = new_names[stage]

        if i == len(new_names) - 1: # Just go back to rigid. Stage before rigid as no linked file
            previous_tp_str = '(InitialTransformParametersFileName "NoInitialTransform")\n'

        else:
            previous_tp = f'{list(tforms)[len(new_names) - (i + 2) ]}.txt'
            previous_tp_str ='(InitialTransformParametersFileName "{}")\n'.format(previous_tp)

        with open(tp, 'r') as fh:
            lines = []

            for line in fh:

                if line.startswith('(InitialTransformParametersFileName'):
                    lines.append(previous_tp_str)
                    continue
                lines.append(line)

        with open(tp, 'w') as wh:
            for line in lines:
                wh.write(line)

    return list(new_names.values())[-1]


def _compare_points(target_transformed_points: Path, moving_points: Path) -> pd.DataFrame:
    """
    For each coresponding pairs of points in points_1 and points_2 compute the euclidean distance

    Parameters
    ----------
    target_transformed_points
        path to an elastix format outputpoints file
    moving_points
        path to Slicer format input points file

    Returns
    -------
    DataFrame with distance between corresponding points
    """

    df_mov = Points().set_slicer(moving_points).points
    df_target_trans = Points().set_elastix_output(target_transformed_points).points

    if len(df_mov) != len(df_target_trans):
        raise ValueError('Fixed and moving should have the same nuber of points')

    results = []
    # Just until I work out how to reshape the input arrays for pdist correctly
    for i in range(len(df_mov)):
        m = df_mov.iloc[i].abs()
        t = df_target_trans.iloc[i].abs()
        dist = spatial.distance.pdist([m, t])
        results.append(dist)
    df = pd.DataFrame.from_records(results)
    df.columns = ['distance']
    return df


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser("Points transformation and distane calculation")
    parser.add_argument('-d', '--lama_output_dir', dest='lama_dir',
                        help='A Lama output directory for a specimen which must at least contain conatin an "output/registration" dir',
                        required=True, type=Path)
    parser.add_argument('-f', '--fixed_points', dest='fixed_points',
                        help='Slicer fiducial file containing markups on the fixed image', required=True, type=Path)
    parser.add_argument('-m', '--moving_points', dest='moving_points',
                        help='Slicer fiducial file containing markups on the moving image', required=True, type=Path)
    parser.add_argument('-o', '--out_dir', dest='out_dir',
                        help='Optional output directory. If not set will make and put in "output/inverted_points"',
                        required=False, type=Path, default=None)

    args = parser.parse_args()
    run(args.lama_dir, args.fixed_points, args.moving_points, args.out_dir)
