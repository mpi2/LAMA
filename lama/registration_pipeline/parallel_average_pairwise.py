"""
Create populaiton avegae from pairwise regitrations and distribute the jobs.

"""

from pathlib import Path
import traceback
from typing import List
from itertools import permutations
import time
from lama.registration_pipeline.validate_config import LamaConfig
from lama.registration_pipeline.run_lama import generate_elx_parameters, ELX_PARAM_PREFIX
from lama import common
from lama.elastix.elastix_registration import TargetBasedRegistration, PairwiseBasedRegistration
from logzero import logger as logging
import logzero
import SimpleITK as sitk



def mean_transform(pair_reg_dir: Path,
                   pre_avg_dir: Path,
                   current_avg_dir) :

    """
    Apply the mean transform and create a populaiton avergae for each stage (latter for vis purpises only)

    Parameters
    ----------
    moving_vol_dir
        A stage root directory. Subdirectories for each fixed image

    Returns
    -------

    """

    # This loop gets all the avergaes


    fixed_vol = list(pre_avg_dir.glob(f'**/{pair_reg_dir.name}.nrrd'))[0]
    tp_paths = list(pair_reg_dir.glob('**/TransformParameters.0.txt'))

    out_dir = current_avg_dir / pair_reg_dir.name
    out_dir.mkdir(exist_ok=True)

    tp_out_name = f'{pair_reg_dir.name}.txt'

    PairwiseBasedRegistration.generate_mean_tranform(tp_paths, fixed_vol, out_dir, tp_out_name, filetype='nrrd')


def get_pairs(inputs_dir):
    """
    Given the input directory return the pairwise combinations (in the form name1_name_2 minus extensions: (ids tuple))

    """
    specimen_ids = [Path(x).stem for x in common.get_file_paths(inputs_dir)]
    perms = list(permutations(specimen_ids, r=2))
    d = {}
    for p in perms:
        k = f'{p[0]}_{p[1]}'
        d[k] = p
    # k = {(f'{x[0]}_{x[1]}'): x for x in perms}
    return d


def get_next_pair(stage_dir: Path):
    """
    Get the next pair to register. If all started return None
    Parameters
    ----------
    stage_dir

    Returns
    -------

    """

    started_dirs = []
    for fixed_root in stage_dir.iterdir():
        if not fixed_root.is_dir():
            continue
        for pair_dir in fixed_root.iterdir():
            started_dirs.append(pair_dir.name)

    return started_dirs


def do_reg(moving, fixed, stage_dir, pair_name, config, elxparam_path):

    # For each volume keep all the regisrtaitons where it is fixed in this folder
    fixed_out_root = stage_dir / fixed.stem

    pair_dir = fixed_out_root / pair_name
    pair_dir.mkdir(exist_ok=True, parents=True)

    fixed_mask = None

    # Do the registrations
    registrator = TargetBasedRegistration(elxparam_path,
                                          fixed,
                                          pair_dir,
                                          config['filetype'],
                                          config['threads'],
                                          fixed_mask
                                          )

    registrator.set_target(moving)
    registrator.rename_output = False

    registrator.run()


def do_stage_reg(pairs: List[str],
                 stage_status_dir: Path,
                 reg_stage_dir: Path,
                 previous_mean_dir,
                 elastix_param_path,
                 config,
                 first=True):
    """
    Process a stage. Only return when fully complelted

    Returns
    -------

    """
    finsihed_stage = False

    while True:  # Pick up unstarted pairwise registrations. Only break when reg and average complete

        # Check if any specimens left (It's possible the avg is being made but all specimens are registered)
        start_dir = stage_status_dir / 'started'
        finish_dir = stage_status_dir / 'finished'
        failed_dir = stage_status_dir / 'failed'


        started = list([x.name for x in start_dir.iterdir()])
        pair_keys = list(pairs.keys())
        not_started = set(pair_keys).difference(started)

        if len(not_started) > 0:
            next_pair_name = list(not_started)[0]  # Some specimens left. Pick up spec_id and process
            next_pair = pairs[next_pair_name]
            # drop a file in the started dir
            start_file = start_dir / next_pair_name

            try:
                with open(start_file, 'x') as fh:
                    if first:  # The first stage all inputs in same dir
                        moving = previous_mean_dir / f'{next_pair[0]}.nrrd'
                        fixed = previous_mean_dir / f'{next_pair[1]}.nrrd'
                    else:  # Outputs are in individual folders
                        moving = previous_mean_dir / next_pair[0] / f'{next_pair[0]}.nrrd'
                        fixed = previous_mean_dir / next_pair[1] / f'{next_pair[1]}.nrrd'

                    try:
                        do_reg(moving, fixed, reg_stage_dir, next_pair_name, config, elastix_param_path)
                        reg_finished_file = finish_dir / next_pair_name
                        open(reg_finished_file, 'x').close()
                    except Exception as e:
                        failed_reg_file = failed_dir / next_pair_name

                        with open(failed_reg_file, 'x') as fh:
                            traceback.print_exc(file=fh)
                        raise

            except FileExistsError:
                # This is raised if another instance has already touched the file.
                continue

        else:  # All specimens have at least been started

            # This block waits for all registrations to be finished or check for any failaures
            while True:
                finished = [x.name for x in finish_dir.iterdir()]
                failed = [x.name for x in failed_dir.iterdir()]

                if failed:
                    raise RuntimeError('Exiting as a failure has been detected')

                not_finished = set(pairs.keys()).difference(finished)

                if not_finished:
                    print("Wating for specimens to finish reg")
                    time.sleep(5)
                else:
                    finsihed_stage = True
                    break
        if finsihed_stage:
            break


def do_mean_transforms(pairs, stage_status_dir, reg_stage_dir, mean_dir, previous_mean_dir, avg_out):

    mean_started_dir = stage_status_dir / 'mean_started'
    mean_finished_dir = stage_status_dir / 'mean_finished'

    moving_ids = []
    for moving_vol_dir in reg_stage_dir.iterdir():
        if not moving_vol_dir.is_dir():
            continue
        moving_ids.append(moving_vol_dir.name)

        try:
            with open(mean_started_dir / moving_vol_dir.name, 'x'):
                mean_transform(moving_vol_dir, previous_mean_dir, mean_dir)
                mean_finished_file = mean_finished_dir / moving_vol_dir.name
                open(mean_finished_file, 'x').close()
        except FileExistsError:
            continue

    while True:
    # Wait for mean transformas to finish
        means_finshed = [x.name for x in mean_finished_dir.iterdir()]
        means_not_finished = set(moving_ids).difference(means_finshed)

        if len(means_not_finished) > 0:
            print('waiting for mean transfroms')
            time.sleep(2)
        else:
            break
    # make averge images
    avg_started_file = str(avg_out) + 'started'
    try:
        with open(avg_started_file, 'x'):
            img_paths = common.get_file_paths(mean_dir)
            avg = common.average(img_paths)
            sitk.WriteImage(avg, str(avg_out))
    except FileExistsError:
        return # We don't have to wait for avergae to finish as it's not required for the next stage


def job_runner(config_path: Path) -> Path:
    """
    Run the registrations specified in the config file

    Returns
    -------
    The path to the final registrered images
    """

    # Load the lama config
    config = LamaConfig(config_path)

    mean_transforms = config_path.parent / 'mean_transforms'
    mean_transforms.mkdir(exist_ok=True)

    avg_dir = config_path.parent / 'averages'
    avg_dir.mkdir(exist_ok=True)

    # Folder to create logic control files
    status_dir = config_path.parent / 'status'
    status_dir.mkdir(exist_ok=True)

    # Get list of specimens
    inputs_dir = config.options['inputs']
    pairs = get_pairs(inputs_dir)

    previous_mean_dir = inputs_dir
    first = True

    for i, reg_stage in enumerate(config['registration_stage_params']):

        stage_id = reg_stage['stage_id']
        avg_out = avg_dir / f'{stage_id}.nrrd'
        logging.info(stage_id)
        reg_stage_dir = Path(config.stage_dirs[stage_id])

        # Make stage dir if not made by another instance of the script
        reg_stage_dir.mkdir(exist_ok=True, parents=True)

        elastix_stage_parameters = generate_elx_parameters(config, do_pairwise=True)[stage_id]
        # Make the elastix parameter file for this stage
        elxparam_path = reg_stage_dir / f'{ELX_PARAM_PREFIX}{stage_id}.txt'

        if not elxparam_path.is_file():
            with open(elxparam_path, 'w') as fh:
                fh.write(elastix_stage_parameters)

        stage_mean_dir = mean_transforms / stage_id
        stage_mean_dir.mkdir(exist_ok=True, parents=True)

        stage_status_dir = status_dir / stage_id
        stage_status_started = stage_status_dir / 'started'
        stage_status_failed = stage_status_dir / 'failed'
        stage_status_finished = stage_status_dir / 'finished'

        stage_status_started.mkdir(exist_ok=True, parents=True)
        stage_status_failed.mkdir(exist_ok=True, parents=True)
        stage_status_finished.mkdir(exist_ok=True, parents=True)

        stage_tform_started = stage_status_dir / 'mean_started'
        stage_tform_finished = stage_status_dir / 'mean_finished'
        stage_tform_started.mkdir(exist_ok=True)
        stage_tform_finished.mkdir(exist_ok=True)

        do_stage_reg(pairs, stage_status_dir, reg_stage_dir, previous_mean_dir,
                     elxparam_path, config, first)

        do_mean_transforms(pairs, stage_status_dir, reg_stage_dir, stage_mean_dir, previous_mean_dir, avg_out)

        first = False
        previous_mean_dir = stage_mean_dir


if __name__ == '__main__':
    import sys
    config_path_ = Path(sys.argv[1])

job_runner(config_path_)
