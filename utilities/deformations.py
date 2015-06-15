#!/usr/bin/env python


import sys
import os
import subprocess
import shutil
import harwellimglib as hil

def run(indir,outdir):
    jacdir = os.path.join(outdir, 'spatialJacobians')
    defdir = os.path.join(outdir, 'deformationFields')
    if not os.path.isdir(jacdir):
        os.mkdir(jacdir)
    if not os.path.isdir(defdir):
        os.mkdir(defdir)
    generate_deformation_fields(indir, defdir, jacdir, 'nrrd')




def generate_deformation_fields(registration_dir, deformation_dir, jacobian_dir, filetype):
    """
    Run transformix on the specified registration stage to generate deformation fields and spatial jacobians
    :param registration_dir:
    :param deformation_dir:
    :param jacobian_dir:
    :return:
    """
    print('### Generating deformation files ###')
    print registration_dir, deformation_dir, jacobian_dir
    # Get the transform parameters
    dirs_ = os.listdir(registration_dir)
    print dirs_
    for dir_ in dirs_:
        if os.path.isdir(os.path.join(registration_dir, dir_)):

            elx_tform_file = os.path.join(registration_dir, dir_, 'TransformParameters.0.txt')
            if not os.path.isfile(elx_tform_file):
                sys.exit("### Error. Cannot find elastix transform parameter file: {}".format(elx_tform_file))

            cmd = ['transformix',
                   '-tp', elx_tform_file,
                   '-out', deformation_dir,
                   '-def', 'all',
                   '-jac', 'all']

            try:
                output = subprocess.check_output(cmd)
            except subprocess.CalledProcessError as e:
                sys.exit('### Transformix failed ###\nError message: {}\nelastix command:{}'.format(cmd))
            else:
                # Rename the outputs and move to correct dirs
                volid = os.path.basename(dir_)
                deformation_out = os.path.join(deformation_dir, 'deformationField.{}'.format(filetype))
                jacobian_out = os.path.join(deformation_dir, 'spatialJacobian.{}'.format(filetype))
                deformation_renamed = os.path.join(deformation_dir, '{}_deformationFied.{}'.format(volid, filetype))
                jacobian_renamed = os.path.join(jacobian_dir, '{}_spatialJacobian.{}'.format(volid, filetype))
                shutil.move(deformation_out, deformation_renamed)
                shutil.move(jacobian_out, jacobian_renamed)

if __name__ == '__main__':
    indir = sys.argv[1]
    outdir = sys.argv[2]

    run(indir, outdir)