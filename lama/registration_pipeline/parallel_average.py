"""
020620

A go at making populaiton average construction parallelizable across the grid.

Take a root directory with inputs. Poll this directory to get a list of specimens.
Look in the output directory and see what folders are avaialbe to see if there are remaiing specimens to process for a
stage.

If one job runner finds a stage complete, but there's not completion file, create an exlusive file and then go about
starting the next stage and.


example
-------
parallel_average.py lama_config.toml

"""

from pathlib import Path
from os.path import join
import time
from lama.registration_pipeline.validate_config import LamaConfig
from lama.registration_pipeline.run_lama import generate_elx_parameters, ELX_PARAM_PREFIX
from lama import common
from lama.elastix.elastix_registration import TargetBasedRegistration
from logzero import logger as logging
import logzero
import SimpleITK as sitk


def make_avg(root_dir: Path,  out_path: Path, log_path):

    paths = []
    for spec_dir in root_dir.iterdir():
        if not spec_dir.is_dir():
            continue
        paths.append(spec_dir / f'{spec_dir.name}.nrrd')

    avg = common.average(paths)
    logzero.logfile(log_path)
    logging.info(f'\nCreating average from:\n')
    logging.info('\n'.join([str(x) for x in paths]))
    sitk.WriteImage(avg, str(out_path))


def check_stage_done(root_dir) -> bool:
    # Check to see if all specimens have finished within a stage
    for spec_dir in root_dir.iterdir():
        if not spec_dir.is_dir():
            continue
        if not (spec_dir / 'spec_done').is_file():
            return False
    return True


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

    # Set the fixed volume up for the first stage. This will checnge each stage if doing population average
    fixed_vol = config['fixed_volume']

    # Get list of specimens
    inputs_dir = config.options['inputs']
    spec_ids = [Path(x).stem for x in common.get_file_paths(inputs_dir)]

    for i, reg_stage in enumerate(config['registration_stage_params']):

        stage_id = reg_stage['stage_id']
        logging.info(stage_id)
        stage_dir = Path(config.stage_dirs[stage_id])

        # Make stage dir if not made by another instance of the script
        stage_dir.mkdir(exist_ok=True, parents=True)

        starting_avg = stage_dir / 'avg_started'
        average_done = stage_dir / "avg_done"

        while True:  # Pick up unstarted speciemens. Only break when reg and avergae complete

            # Check if any specimens left (It's possible the avg is being made but all specimens are registered)
            spec_stage_dirs = [x.name for x in stage_dir.iterdir() if x.is_dir()]
            not_started = set(spec_ids).difference(spec_stage_dirs)

            next_stage = False  # No breaking out yet

            if len(not_started) > 0:
                next_spec_id = list(not_started)[0] # Some specimens left. Pick up spec_id and process

            else:  # All specimens are being processed
                next_stage = True

                #  This block controls what happens if we have all speciemns registered
                while True:
                    if not check_stage_done(stage_dir):
                        print('waiting for stage to finish')
                        time.sleep(5)
                        continue

                    print('stage finished')

                    if average_done.is_file():
                        print('found average done file')
                        break  # Next stage
                    else:
                        if starting_avg.is_file():
                            print('found starting average file')
                            time.sleep(5)
                            continue
                        else:
                            try:
                                open(starting_avg, 'x')
                            except FileExistsError:
                                time.sleep(5)
                                print('cannot write avg starting file')
                                continue
                            else:
                                average_path = avg_dir / f'{stage_id}.nrrd'
                                make_avg(stage_dir, average_path, avg_dir / f'{stage_id}.log')
                                open(average_done, 'x').close()
                                print('making average')
                                break

            if next_stage:
                print('breaking stage')
                break

            # Get the input for this specimen
            if i == 0:  # The first stage
                moving = inputs_dir / f'{next_spec_id}.nrrd'
            else:
                moving = list(config.stage_dirs.values())[i-1] / next_spec_id / f'{next_spec_id}.nrrd'
                fixed_vol = avg_dir / f'{list(config.stage_dirs.keys())[i-1]}.nrrd'
            reg_method = TargetBasedRegistration

            # Make the elastix parameter file for this stage
            elxparam = elastix_stage_parameters[stage_id]
            elxparam_path = stage_dir / f'{ELX_PARAM_PREFIX}{stage_id}.txt'

            if not elxparam_path.is_file():
                with open(elxparam_path, 'w') as fh:
                    if elxparam:
                        fh.write(elxparam)

            fixed_mask = None

            logging.info(moving)

            # Do the registrations
            registrator = reg_method(elxparam_path,
                                     moving,
                                     stage_dir,
                                     config['filetype'],
                                     config['threads'],
                                     fixed_mask
                                     )

            registrator.set_target(fixed_vol)

            try:
                registrator.run()  # Do the registrations for a single stage
            except FileExistsError as e:
                # 040620: Bodge as some specimens are picked up twice.
                # Need a better way to make sure each speciemn picked up only once
                continue

            spec_done = stage_dir / next_spec_id / 'spec_done'  # The directory gets created in .run()
            open(spec_done, 'x').close()


if __name__ == '__main__':
    import sys
    config_path_ = Path(sys.argv[1])

job_runner(config_path_)