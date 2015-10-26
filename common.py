import logging
import subprocess
import os
from os.path import join, basename
import sys
import datetime
import shutil
import SimpleITK as sitk

LOG_MODE = logging.DEBUG


def write_array(array, path):
    """
    Write a numpy array to and image file using SimpleITK
    """
    sitk.WriteImage(sitk.GetImageFromArray(array), path)


def img_path_to_array(img_path):
    if os.path.isfile(img_path):
        try:
            img = sitk.ReadImage(img_path)

        except RuntimeError:
            print "Simple ITK cannot read {}".format(img_path)
            return None
        else:
            array = sitk.GetArrayFromImage(img)
            return array
    else:
        print '{} is not a real path'.format(img_path)
        return None

def init_log(logpath, name):

    logging.basicConfig(filename=logpath, level=LOG_MODE, filemode="w")
    logging.info(name)

    # switch working dir to this module's location
    orig_wd = os.getcwd()
    module_dir = os.path.dirname(os.path.realpath(__file__))
    os.chdir(module_dir)

    try:
        git_log = subprocess.check_output(['git', 'log', '-n', '1'])
        git_commit = git_log.splitlines()[0]
        git_branch = subprocess.check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD'])
    except subprocess.CalledProcessError:
        git_commit = "Git commit info not available"
        git_branch = "Git branch not available"

    logging.info("git branch: {}".format(git_branch.strip()))
    logging.info("git {}".format(git_commit.strip()))

    # back to original working dir
    os.chdir(orig_wd)

def log_time(msg):
    now = datetime.datetime.now()
    logging.info("{}: {}/{}/{} - {}:{}".format(msg, now.day, now.month, now.year, now.hour, now.minute))

def mkdir_force(dir_):
    if os.path.isdir(dir_):
        shutil.rmtree(dir_)
    os.mkdir(dir_)

def mkdir_if_not_exists(dir_):
    if not os.path.exists(dir_):
        os.makedirs(dir_)


def GetFilePaths(folder, extension_tuple=('.nrrd', '.tiff', '.tif', '.nii', '.bmp', 'jpg', 'mnc', 'vtk'), pattern=None):
    """
    Test whether input is a folder or a file. If a file or list, return it.
    If a dir, return all images within that directory.
    Optionally test for a pattern to sarch for in the filenames
    """
    if not os.path.isdir(folder):
        if isinstance(folder, basestring):
            return [folder]
        else:
            return folder
    else:
        paths = []
        for root, _, files in os.walk(folder):
            for filename in files:
                if filename.lower().endswith(extension_tuple):
                    if pattern:
                        if pattern and pattern not in filename:
                            continue
                    #paths.append(os.path.abspath(os.path.join(root, filename))) #Broken on shared drive
                    paths.append(os.path.join(root, filename))
        return paths


def check_config_entry_path(dict_, key):
    try:
        value = dict_[key]
    except KeyError:
        raise Exception("'{} not specified in config file".format(key))
    else:
        if not os.path.isdir(value):
            raise OSError("{} is not a correct directory".format(value))

