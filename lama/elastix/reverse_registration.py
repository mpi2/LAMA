"""
This module provides the same function as lama.elastix.invert_transforms.py in that it creates Transform parameters
that can be used to propagate labels from atlas to moving image.

The diffrence is in the approach. In this module istead of minimizing the DiplacementmagnitudePenalty using the forawrd
transform as shown in the elastix manual, we do a registration of the population avg -> image of iterest than use the
transforms genrted to propagate the labels.


Example structure of output folder
-----------------------------------

.
├── invert.yaml
├── affine
│   └── 1574928_download
│       ├── elastix.log
│       ├── ImageInvertedTransform.txt
│       ├── inversion_parameters.txt
│       ├── IterationInfo.0.R0.txt
│       ├── IterationInfo.0.R1.txt
│       ├── labelInvertedTransform.txt
│       └── TransformParameters.0.txt
├── deformable_128
│   └── 1574928_download
│       ├── elastix.log
│       ├── ImageInvertedTransform.txt
│       ├── inversion_parameters.txt
│       ├── IterationInfo.0.R0.txt
│       ├── labelInvertedTransform.txt
│       └── TransformParameters.0.txt
└── rigid
    └── 1574928_download
        ├── elastix.log
        ├── ImageInvertedTransform.txt
        ├── inversion_parameters.txt
        ├── IterationInfo.0.R0.txt
        ├── IterationInfo.0.R1.txt
        ├── labelInvertedTransform.txt
        └── TransformParameters.0.txt



"""

from pathlib import Path
import tempfile
import os
import shutil
import subprocess
from collections import defaultdict
from multiprocessing import Pool
from os.path import join, abspath, isfile
from typing import Union, List, Dict
import os

from logzero import logger as logging
import yaml
import yaml

from lama import common
from lama.common import cfg_load
from lama.registration_pipeline.validate_config import LamaConfig
from lama.registration_pipeline import run_lama
from lama.elastix.invert_transforms import get_reg_dirs

from lama.elastix import (ELX_TRANSFORM_NAME, ELX_PARAM_PREFIX, LABEL_INVERTED_TRANFORM,
                          IMAGE_INVERTED_TRANSFORM, INVERT_CONFIG, RESOLUTION_IMGS_DIR, IMG_PYRAMID_DIR)
from lama.elastix.invert_transforms import (LABEL_REPLACEMENTS, IMAGE_REPLACEMENTS,
                                            )
from lama.elastix.elastix_registration import TargetBasedRegistration


def reverse_registration(config: Union[str, LamaConfig]):
    common.test_installation('elastix')

    if isinstance(config, (Path, str)):
        config = LamaConfig(Path(config))

    # threads = str(config['threads'])

    inv_outdir = config.mkdir('inverted_transforms')

    # Set the moving volume dir and the fixed image for the first stage
    fixed_vols_dir = config['inputs']
    # Get the fixed vols
    fixed_vol_paths = common.get_file_paths(fixed_vols_dir)

    # Get the fixed and moving images. They are flipped compared to the forward registration
    moving_vol = config['fixed_volume']

    # Do the forward registration for each imaage (usually just one image if using the jobrunner script)
    for fixed_vol in fixed_vol_paths:
        run_registration_schedule(config, fixed_vol, moving_vol, inv_outdir)


def run_registration_schedule(config: LamaConfig, fixed_vol, moving_vol: Path, outdir: Path) -> Path:
    """
    Run the registrations specified in the config file, but flip the moving and fixed images

    Returns
    -------
    The path to the final registrered images
    """

    elastix_stage_parameters = run_lama.generate_elx_parameters(config)

    # Set the moving volume dir and the fixed image for the first stage
    #  Set the fixed volume up for the first stage. This will checnge each stage if doing population average
    stage_ids = []
    for i, reg_stage in enumerate(config['registration_stage_params']):

        tform_type = reg_stage['elastix_parameters']['Transform']
        euler_stage = True if tform_type == 'EulerTransform' else False

        #  Make the stage output dir
        stage_id = reg_stage['stage_id']
        stage_ids.append(stage_id)
        stage_dir = outdir / stage_id
        stage_dir.mkdir(exist_ok=True)

        logging.info("### Reverse registration - current registration step: {} ###".format(stage_id))

        # Make the elastix parameter file for this stage
        elxparam = elastix_stage_parameters[stage_id]
        elxparam_path = join(stage_dir, ELX_PARAM_PREFIX + stage_id + '.txt')

        with open(elxparam_path, 'w') as fh:
            fh.write(elxparam)

        # maybe we should add fixed mask
        fixed_mask = None

        # Do the registrations
        registrator = TargetBasedRegistration(elxparam_path,
                                 moving_vol,
                                 stage_dir,
                                 config['filetype'],
                                 config['threads'],
                                 fixed_mask
                                 )

        registrator.set_target(fixed_vol)



        if reg_stage['elastix_parameters']['Transform'] == 'BSplineTransform':
            logging.info(f'Folding correction for stage {stage_id} set')
            registrator.fix_folding = config['fix_folding']  # Curently only works for TargetBasedRegistration

        registrator.run()  # Do the registrations for a single stage
        os.remove(elxparam_path)

        # As the stage output diretory is named as the moving image, but ion the case we eant it named the same as the
        # fixed image
        stage_spec_dir = next(stage_dir.glob(f'*{moving_vol.stem}'))
        new_stage_spec_dir = stage_dir / fixed_vol.stem
        stage_spec_dir.rename(new_stage_spec_dir)

        # Now delete everything we don't need
        to_keep = [ELX_TRANSFORM_NAME, 'elastix.log']
        for f in new_stage_spec_dir.iterdir():
            if f.name not in to_keep:
                try:
                    shutil.rmtree(f)
                except NotADirectoryError:
                    f.unlink()

        src_tform_file = stage_dir / fixed_vol.stem / ELX_TRANSFORM_NAME
        label_tform_file = stage_dir / fixed_vol.stem / LABEL_INVERTED_TRANFORM
        image_tform_file = stage_dir / fixed_vol.stem / IMAGE_INVERTED_TRANSFORM
        modify_elx_parameter_file(src_tform_file, label_tform_file, LABEL_REPLACEMENTS)
        modify_elx_parameter_file(src_tform_file, image_tform_file, IMAGE_REPLACEMENTS)



    logging.info("### Reverse registration finished ###")

    # Create invert.yaml. For reverse registration based inversion, it will be the same order as the registrations
    # were done
    d = {'inversion_order': stage_ids}
    with open(outdir / 'invert.yaml', 'w') as fh:
        yaml.dump(d, fh)





def modify_elx_parameter_file(elx_param_file: Path, newfile_name: str, replacements: Dict):
    """
    Modifies the elastix input parameter file that was used in the original transformation.
    Adds DisplacementMagnitudePenalty (which is needed for inverting)
    Turns off writing the image results at the end as we only need an inverted output file.
    Also changes interpolation order in the case of inverting labels

    Parameters
    ----------
    elx_param_file
        path to elastix input parameter file
    newfile_name
        path to save modified parameter file to
    replacements


    """

    try:
        with open(elx_param_file) as old, open(newfile_name, "w") as new:

            for line in old:
                if line.startswith('(WriteResultImage '):
                    line = '(WriteResultImage "false")\n'
                if line.startswith('WriteResultImageAfterEachResolution '):
                   continue
                try:
                    param_name = line.split()[0][1:]
                except IndexError:
                    continue  # comment?

                if param_name in replacements:
                    value = replacements[param_name]
                    try:
                        int(value)
                    except ValueError:
                        # Not an int, neeed quotes
                        line = '({} "{}")\n'.format(param_name, value)
                    else:
                        # An int, no quotes
                        line = '({} {})\n'.format(param_name, value)
                new.write(line)
    except IOError as e:
        logging.error("Error modifying the elastix parameter file: {}".format(e))
        return False
    return True





if __name__ == '__main__':
    import sys
    cfg_path = sys.argv[1]
    reverse_registration(cfg_path)