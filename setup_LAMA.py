#! /usr/bin/env python3

import subprocess as sub
import os

# How to ensure we use python3.6?

VENV_NAME = 'lama-venv'
# Intall packages via pip3

# Add the .pth file into sitepackages so we can load in modules

script_dir = os.path.dirname(os.path.abspath(__file__))

venv_path = os.path.join(script_dir,  VENV_NAME)

print('installing python virtual environment')

try:
    res = sub.check_output(['python3', '-m', 'venv', venv_path])
except sub.CalledProcessError as e: # Migt
    print(e.output)
    raise SystemExit
python_bin = os.path.join(script_dir, VENV_NAME, 'bin', 'python3')
pip3_bin = os.path.join(script_dir, VENV_NAME, 'bin', 'pip3')

req = os.path.join(script_dir, 'requirements.txt')

print('installing python packages')
sub.check_output([python_bin, pip3_bin, 'install', '-r', req])

site_packages = os.path.join(script_dir, VENV_NAME, 'lib', 'python3.6', 'site-packages')

pth_file_path = os.path.join(site_packages, 'lama.pth')

dirs_to_add_to_path = ['lib', 'elastix', 'img_processing']

print('Setting up paths')
with open(pth_file_path, 'w') as fh:
    for dir_ in dirs_to_add_to_path:
        path_to_add = os.path.join(script_dir, dir_)
        fh.write("{}\n".format(path_to_add))

print('Succesfully setup virtual envoronment\n'
      "run 'source lama-venv/bin/activate'\n"
      "Then to run lama\npython3 ./run_lama.py")
