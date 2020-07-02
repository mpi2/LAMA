"""
Create populaiton avegae from pairwise regitrations and distribute the jobs.

Note:
"""

from pathlib import Path
from os.path import join
from itertools import permutations
import time
from lama.registration_pipeline.validate_config import LamaConfig
from lama.registration_pipeline.run_lama import generate_elx_parameters, ELX_PARAM_PREFIX
from lama import common
from lama.elastix.elastix_registration import TargetBasedRegistration, PairwiseBasedRegistration
from logzero import logger as logging
import logzero
import SimpleITK as sitk
import subprocess as sub


def check_stage_done(root_dir) -> bool:
    # Check to see if all specimens have finished within a stage
    for spec_dir in root_dir.iterdir():
        if not spec_dir.is_dir():
            continue
        if not (spec_dir / 'spec_done').is_file():
            return False
        if (spec_dir / 'failed').is_file():
            raise RuntimeError(f'Reg failure: {spec_dir}')
    return True


def transorm(root_reg_dir: Path,
             pre_avg_dir: Path,
             current_avg_dir) :

    """
    Apply the mean transform and cerate a populaiton avergae for each stage (latter for vis purpises only)

    Parameters
    ----------
    root_reg_dir
        A stage root directory. Subdirectories for each fixed image

    Returns
    -------

    """
    paths = []

    # This loop gets all the avergaes
    for fixed_dir in root_reg_dir.iterdir():

        if not fixed_dir.is_dir():
            continue

        # Check that all reg is finished
        while True:
            if check_stage_done(fixed_dir):
                break
            else:
                print('Waitng for reg to finish')
                time.sleep(5)

        avg_started = fixed_dir / 'avg_started'
        avg_finished = fixed_dir / 'finished'

        try:
            open(avg_started, 'x').close()
        except FileExistsError:
            continue

        fixed_vol = list(pre_avg_dir.glob(f'**/{fixed_dir.name}.nrrd'))[0]
        tp_paths = list(fixed_dir.glob('**/TransformParameters.0.txt'))

        out_dir = current_avg_dir / fixed_dir.name
        out_dir.mkdir(exist_ok=True)

        tp_out_name = f'{fixed_dir.name}.txt'

        PairwiseBasedRegistration.generate_mean_tranform(tp_paths, fixed_vol, out_dir, tp_out_name, filetype='nrrd')
        open(avg_finished, 'x').close()

    # At this point we have looped over each fixed root dir and either setup avergae creation or skippd if
    # already started. Now wait until all completed and return
    for fixed_dir in root_reg_dir.iterdir():#
        if not fixed_dir.is_dir():
            continue
        while True:
            if avg_finished.is_file():
            # if (fixed_dir / 'finished').is_file():
                break
            else:
                print(f'waiting for {avg_finished}')
                time.sleep(5)
    print('finshed transforms')
    return True






def get_pairs(inputs_dir):
    """
    Given the input directory return the pairwise combunations (in the form name1_name_2 minus extensions: (ids tuple))

    """
    specimen_ids = [Path(x).stem for x in common.get_file_paths(inputs_dir)]
    perms = list(permutations(specimen_ids, r=2))
    k = {(f'{x[0]}_{x[1]}'): x for x in perms}
    return k


def get_started_pairwse_reg(stage_dir: Path):

    started_dirs = []
    for fixed_root in stage_dir.iterdir():
        if not fixed_root.is_dir():
            continue
        for pair_dir in fixed_root.iterdir():
            started_dirs.append(pair_dir.name)

    return started_dirs


def job_runner(config_path: Path) -> Path:
    """
    Run the registrations specified in the config file

    Returns
    -------
    The path to the final registrered images
    """

    config = LamaConfig(config_path)
    print(common.git_log())


    avg_dir = config.options['average_folder']
    avg_dir.mkdir(exist_ok=True, parents=True)

    elastix_stage_parameters = generate_elx_parameters(config, do_pairwise=config['pairwise_registration'])

    # Get list of specimens
    inputs_dir = config.options['inputs']
    pairs = get_pairs(inputs_dir)

    previous_av_dir = inputs_dir

    for i, reg_stage in enumerate(config['registration_stage_params']):

        stage_id = reg_stage['stage_id']
        logging.info(stage_id)
        stage_dir = Path(config.stage_dirs[stage_id])

        # Make stage dir if not made by another instance of the script
        stage_dir.mkdir(exist_ok=True, parents=True)

        stage_avg_dir = avg_dir / stage_id
        stage_avg_dir.mkdir(exist_ok=True, parents=True)

        while True:  # Pick up unstarted pairwise registrations. Only break when reg and average complete

            # Check if any specimens left (It's possible the avg is being made but all specimens are registered)
            pairwisre_stage_dirs = get_started_pairwse_reg(stage_dir)
            not_started = set(pairs).difference(pairwisre_stage_dirs)

            next_stage = False  # No breaking out yet

            if len(not_started) > 0:
                next_pair_name = list(not_started)[0] # Some specimens left. Pick up spec_id and process
                next_pair = pairs[next_pair_name]

            else:  # All specimens are being processed
                next_stage = True

                #  This block controls what happens if we have all speciemns registered
                while True:
                    if not transorm(stage_dir, previous_av_dir, stage_avg_dir): # returns True when all averages are finished
                        print('waiting for stage to finiish')
                        time.sleep(5)
                        continue
                    else:
                        break

            if next_stage:
                print('breaking stage')
                previous_av_dir = stage_avg_dir
                break

            # Get the input for this specimen
            if i == 0:  # The first stage
                moving = inputs_dir / f'{next_pair[0]}.nrrd'
                fixed = inputs_dir /  f'{next_pair[1]}.nrrd'
            else:
                moving = previous_av_dir / next_pair[0] / f'{next_pair[0]}.nrrd'
                fixed = previous_av_dir / next_pair[1] / f'{next_pair[1]}.nrrd'

            # Make the elastix parameter file for this stage
            elxparam = elastix_stage_parameters[stage_id]
            elxparam_path = stage_dir / f'{ELX_PARAM_PREFIX}{stage_id}.txt'

            # For each volume keep all the regisrtaitons where it is fixed in this folder
            fixed_out_root = stage_dir / fixed.stem

            pair_dir = fixed_out_root / next_pair_name
            pair_dir.mkdir(exist_ok=True, parents=True)

            if not elxparam_path.is_file():
                with open(elxparam_path, 'w') as fh:
                    if elxparam:
                        fh.write(elxparam)

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

            try:
                spec_started = pair_dir / 'spec_started'
                open(spec_started, 'x').close()
                registrator.run()  # Do the registrations for a single stage
            except FileExistsError as e:
                #  make sure each specimen picked up only once
                continue
            except Exception:
                failed = pair_dir / 'failed'
                open(failed, 'x').close()
                raise

            spec_done = pair_dir / 'spec_done'  # The directory gets created in .run()
            open(spec_done, 'x').close()


if __name__ == '__main__':
    import sys
    config_path_ = Path(sys.argv[1])

job_runner(config_path_)

# transorm(Path('/mnt/IMPC_media/LAMA_staging/e15_5/080620_pop_avg/1/290620-good-contrast-specimens/output/registrations/affine'),
#          Path('/mnt/IMPC_media/LAMA_staging/e15_5/080620_pop_avg/1/290620-good-contrast-specimens/rigid'),
#          Path('/mnt/IMPC_media/LAMA_staging/e15_5/080620_pop_avg/1/290620-good-contrast-specimens/output/averages/affine'))