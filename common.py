import logging
import subprocess as sub
import shutil
from traceback import format_exception
import SimpleITK as sitk
import os
import psutil
import sys
import numpy as np
import csv
import yaml
from os.path import abspath, join, basename, splitext
from collections import defaultdict, namedtuple, OrderedDict
import re

INDV_REG_METADATA = 'reg_metadata.yaml'

LOG_FILE = 'LAMA.log'
LOG_MODE = logging.DEBUG
DEFAULT_VOXEL_SIZE = 28.0
IMG_EXTS = ['.nii', '.nrrd', '.tif', '.tiff', '.mhd', '.mnc']


Roi = namedtuple('Roi', 'x1 x2 y1 y2 z1 z2')


STAGING_INFO_FILENAME = 'staging_info.csv'
class RegistrationException(Exception):
    """
    An exception that is raised when the current process (inversion, stats etc cannot complete due to problems with the
    data
    """
    pass

class TransformixException(Exception):
    pass


class LamaDataException(Exception):
    """
    An exception that is raised when the current process (inversion, stats etc cannot complete due to problems with the
    data
    """
    pass


def excepthook_overide(exctype, value, traceback):
    """
    USed to override sys.xcepthook so we can log any ancaught Exceptions
    :param exctype:
    :param value:
    :param tb:
    :return:
    """

    if isinstance(exctype, type(LamaDataException)):
        print(''.join(format_exception(exctype, value, traceback)))
        print("\n\n\n\n")
        print('Lama encountered a problem with reading or interpresting some data. Plese check the log files')

    else:
        print("\n\nLAMA encountered an unknown error\nPlease email us at sig.har.mrc.ac.uk with the contents of the LAMA.log\n")
        print(''.join(format_exception(exctype, value, traceback)))


class LoadImage(object):
    def __init__(self, img_path):
        self.img_path = img_path
        self.error_msg = None
        self.img = None
        self._read()

    def __nonzero__(self):
        """
        Overload this so we can do simple 'is LoadImage' to check if img loaded
        """
        if self.img is None:
            return False
        else:
            return True

    @property
    def array(self):
        return sitk.GetArrayFromImage(self.img)

    def _read(self):
        if os.path.isfile(self.img_path):
            try:
                self.img = sitk.ReadImage(self.img_path)
            except RuntimeError:
                self.error_msg = "possibly corrupted file {}".format(self.img_path)
        else:
            self.error_msg = "path does not exist: {}".format(self.img_path)


class PathToITKImage(object):
    def __init__(self, img_path):
        self.img_path = img_path
        self.error_msg = None
        self._to_array()

    def _to_array(self):
        if os.path.isfile(self.img_path):
            try:
                img = sitk.ReadImage(self.img_path)
            except RuntimeError:
                self.error_msg = "possibly corrupted file {}".format(self.img_path)
                return None
            else:
                return img
        else:
            self.error_msg = "path does not exist: {}".format(self.img_path)
            return None


def plot_swarm():
    pass


def write_array(array, path, compressed=True):
    """
    Write a numpy array to and image file using SimpleITK
    """
    sitk.WriteImage(sitk.GetImageFromArray(array), path, compressed)


def img_path_to_array(img_path):
    if os.path.isfile(img_path):
        try:
            img = sitk.ReadImage(img_path)

        except RuntimeError as e:
            raise OSError(e)
        else:
            array = sitk.GetArrayFromImage(img)
            return array
    else:
        return None


def git_log():

    # switch working dir to this module's location
    orig_wd = os.getcwd()
    module_dir = os.path.dirname(os.path.realpath(__file__))
    os.chdir(module_dir)

    try:
        log = sub.check_output(['git', 'log', '-n', '1'])
        git_commit = log.splitlines()[0]
        git_branch = sub.check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD'])
    except sub.CalledProcessError:
        git_commit = "Git commit info not available"
        git_branch = "Git branch not available"

    message = "git branch: {}. git {}: ".format(git_branch.strip(), git_commit.strip())

    # back to original working dir
    os.chdir(orig_wd)
    return message


def init_logging(logpath):

    if os.path.exists(logpath):  # Create new log file if one already exists
        i = 1
        while True:
            path, ext = os.path.splitext(logpath)
            newname = path + '_' + str(i)
            new_logpath = newname + ext
            if not os.path.exists(new_logpath):
                logpath = new_logpath
                break
            i += 1

    logging.basicConfig(filename=logpath, level=LOG_MODE,
                        format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p')

    stdout_log = logging.StreamHandler(sys.stdout)
    return logging.getLogger().addHandler(stdout_log)


def load_label_map_names(organ_names_path, include_terms=False):
    """
    Lod label name csv

    example file format

        label_name,term
        1,l1,emap:1
        2,l2,emap:2
        3,l3,emap:4
        4,l4,emap:5
        5,l5,emap:6


    returns
    -------
    pandas data frame
        label_num,label_name,emapa_term(optional)

    """
    import pandas as pd
    df = pd.read_csv(organ_names_path)

    return df


def mkdir_force(dir_):
    shutil.rmtree(dir_, ignore_errors=True)
    os.mkdir(dir_)


def mkdir_if_not_exists(dir_):
    if not os.path.exists(dir_):
        os.makedirs(dir_)


def get_file_paths(folder, extension_tuple=('.nrrd', '.tiff', '.tif', '.nii', '.bmp', 'jpg', 'mnc', 'vtk', 'bin'),
                   pattern=None, ignore_folder=""):
    """
    Test whether input is a folder or a file. If a file or list, return it.
    If a dir, return all images within that directory.
    Optionally test for a pattern to sarch for in the filenames
    """

    if not os.path.isdir(folder):
        return False
    else:
        paths = []
        for root, subfolders, files in os.walk(folder):
            if ignore_folder in subfolders:
                subfolders.remove(ignore_folder)
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


def print_memory(fn):
    def wrapper(*args, **kwargs):
        process = psutil.Process(os.getpid())
        start_rss, start_vms = process.get_memory_info()
        try:
            return fn(*args, **kwargs)
        finally:
            end_rss, end_vms = process.get_memory_info()
            print((end_rss - start_rss), (end_vms - start_vms))
    return wrapper


def get_inputs_from_file_list(file_list_path, config_dir):
    """
    Get registration inputs from a file

    example filelist:
        dir: relative_path_to_img_folder
        folder_with_img_1
        folder_with_img_2

    Parameters
    ----------
    file_list_path: str
        path to img file list
    config_dir
        path to the file file name list
    Raises
    ------
    OsError
        if file does not exist

    Returns
    -------

    """
    filtered_paths = []
    if not os.path.isfile(file_list_path):
        return None
    with open(file_list_path, 'r') as reader:
        root_path_dict = defaultdict(list)
        root = None
        for line in reader:
            if line.startswith('dir:'):
                root = abspath(join(config_dir, line.strip('dir:').strip()))
                continue
            if not root:
                raise(LamaDataException('The root directory is missing in the image directory list file {}\n'
                                        'first line should contain "dir:relative/path/to/folder/with/images" '.format(file_list_path)))

            base = line.strip()
            root_path_dict[root].append(base)
            i = 0
    for root, bases in root_path_dict.items():

        # if it's an image path load it. If a directory, load all images from it
        for base in bases:
            i += 1
            path = join(root, base)
            if os.path.isdir(path):
                img_paths = get_file_paths(path)
                filtered_paths.extend([abspath(x) for x in img_paths if splitext(basename(x))[0].strip('seg_') in bases])
            else:
                filtered_paths.append(path)

    return filtered_paths


def Average(img_dirOrList, search_subdirs=True):
    '''
    Create an average volume from multiple volumes
    @return: sitk Image
    '''

    if isinstance(img_dirOrList, basestring):
        images = get_file_paths(img_dirOrList)
    else:
        images = img_dirOrList

    # sum all images together
    summed = sitk.GetArrayFromImage(sitk.ReadImage(images[0]))
    for image in images[1:]:  # Ommit the first as we have that already
        np_array = sitk.GetArrayFromImage(sitk.ReadImage(image))
        try:
            summed += np_array
        except ValueError as e:
            print("Numpy can't average this volume {0}".format(image))

    # Now make average
    summed /= len(images)
    avg_img = sitk.GetImageFromArray(summed)
    return avg_img


def subsample(array, chunk_size, mask=False):
    """
    Parameters
    ----------
    array: numpy.ndarray
    """

    shape = array.shape
    i = 0
    subsampled_array = []
    out_shape = [0, 0, 0]  # zyx
    for z in range(0, shape[0] - chunk_size, chunk_size):
        out_shape[0] += 1
        for y in range(0, shape[1] - chunk_size, chunk_size):
            if z == 0:
                out_shape[1] += 1
            for x in range(0, shape[2] - chunk_size, chunk_size):
                if z == 0 and y == 0:
                    out_shape[2] += 1
                mask_region = array[z: z + chunk_size, y: y + chunk_size, x: x + chunk_size]
                if mask:  # If any in the cube is a mask element, make the whole cube a mask
                    if np.any(mask_region):
                        subsampled_array.insert(i, 1)
                    else:
                        subsampled_array.insert(i, 0)
                else:
                    subsampled_array.insert(i, np.mean(mask_region))
                i += 1
    if mask:
        return np.array(subsampled_array).astype(np.bool).reshape(out_shape)
    else:
        return np.array(subsampled_array).reshape(out_shape)


def write_file_list(root_names_dict, outpath):
    """
    Given a dict of root folder->[list of basenames] create a list that can be read by other componabt of the pipleine
    Parameters
    ----------
    root_names_dict
    """
    with open(outpath) as fh:
        for root, basenames in root_names_dict.iteritems():
            fh.write({'dir:{}\n'.format(root)})
            for base in basenames:
                fh.write('{}\n'.format(base))


def csv_read_lines(path):
    """
    Read lines from a csv
    ----------
    path: str
        path to csv file

    Returns
    -------
    list of lines
    """
    lines = []
    with open(path, 'rb') as fh:
        reader = csv.reader(fh)
        for line in reader:
            lines.append(line[0])
    return lines


def csv_read_dict(path):
    """
    Read lines from a csv
    ----------
    path: str
        path to csv file

    Returns
    -------
    dict where column1 = key and column2 = value
    """
    lines = {}
    with open(path, 'rb') as fh:
        reader = csv.reader(fh)
        for line in reader:
            lines[line[0]] = line[1]
    return lines


def select_subset(paths, subset_ids):
    """
    Trim the files found in the wildtype input directory to thise in the optional subset list file
    """
    wt_paths_to_use = []

    for path in paths:
        vol_name = os.path.splitext(os.path.basename(path))[0]
        if vol_name in subset_ids:
            wt_paths_to_use.append(path)
    return wt_paths_to_use


def check_file_paths(paths, ret_string=False):
    """
    Check for path's existence. Return True if all exist. or a list of failed paths 
    Parameters
    ----------
    paths
    ret_string: boolean
        return failed paths as a string one path on each line
    Returns
    -------

    """
    failed = []
    for path in paths:
        if not os.path.isfile(path):
            failed.append(path)
    if not failed:
        return True
    else:
        if ret_string:
            return "\n".join(failed)


class bcolors:
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def specimen_ids_from_paths(file_names):
    """
    Given a list of file paths get the basename minus the img extension, whicj should be the specimen name
    Parameters
    ----------
    file_names

    Returns
    -------

    """
    return [strip_img_extension(basename(x)) for x in file_names]


def specimen_id_from_file_path(file_path):
    return strip_img_extension(basename(file_path))

def strip_img_extensions(file_names):
    result = []
    for f in file_names:
        result.append(strip_img_extension(f))
    return result


def strip_img_extension(file_name):
    ext = splitext(file_name.lower())[1]
    if ext in IMG_EXTS:
        return file_name.rstrip(ext)
    else:
        return file_name


def test_installation(app):
    try:
        sub.check_output([app])
    except Exception:  # can't seem to log CalledProcessError
        logging.error('It looks like {} may not be installed on your system\n'.format(app))
        return False
    else:
        return True


def is_r_installed():
    installed = True
    FNULL = open(os.devnull, 'w')
    try:
        sub.call(['Rscript'], stdout=FNULL, stderr=sub.STDOUT)
    except sub.CalledProcessError:
        installed = False
    except OSError:
        installed = False
        logging.warn('R or Rscript not installed. Will not be able to use linear model')
    return installed
