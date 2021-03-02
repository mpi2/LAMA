import subprocess as sub
import shutil
from pathlib import Path
from traceback import format_exception
from os.path import abspath, join, basename, splitext
from collections import defaultdict, namedtuple
import sys
import os
from datetime import datetime
from typing import Union, List, Tuple, Dict
import urllib, io
import urllib.request
import http.client as http
import zipfile
import csv
import datetime
import signal
import decimal

from logzero import logger as logging
import logzero
import SimpleITK as sitk
import numpy as np
import pandas as pd
import psutil
import argparse

import yaml
import toml

from lama.img_processing import read_minc
from lama.elastix import RESOLUTION_IMGS_DIR, IMG_PYRAMID_DIR
import lama
INDV_REG_METADATA = 'reg_metadata.yaml'

LOG_FILE = 'LAMA.log'
DEFAULT_VOXEL_SIZE = 28.0
IMG_EXTS = ('.nii', '.nrrd', '.tif', '.tiff', '.mhd', '.mnc', '.npy')
INVERTED_MASK_DIR = 'inverted_stats_masks'


Roi = namedtuple('Roi', 'x1 x2 y1 y2 z1 z2')

ORGAN_VOLUME_CSV_FILE = 'organ_volumes.csv'
STAGING_INFO_FILENAME = 'staging_info_volume.csv'
FOLDING_FILE_NAME = 'folding_report.csv'

lama_root_dir = Path(lama.__file__).parent


class CheckSinglePathAction(argparse.Action):
    """Resolve paths during argument parsing and check for existence"""

    def __call__(self, parser, namespace, values, option_string=None):

        path = Path(values).resolve()

        if not path.exists():
            raise argparse.ArgumentTypeError(f"{path} does not exit")

        setattr(namespace, self.dest, path)



def date_dhm() -> datetime.datetime:
    """
    Get the current date stamp minus seconds and microseconds
    Returns
    -------

    """
    return datetime.datetime.now().replace(second=0, microsecond=0)


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


def disable_warnings_in_docker():
    """
    If running in a Docker container, switch off warnings as we are getting cython warnings coming from sklearn and pandas etc.
    If running in Docker, there will be a lama forlder in the root directory
    """

    if os.path.isdir(os.path.join(os.sep, 'lama')):

        print('Lama running in docker')

        import warnings
        warnings.filterwarnings('ignore')


def excepthook_overide(exctype, value, traceback):
    """
    Used to override sys.excepthook so we can log any ancaught Exceptions

    Parameters
    ----------
    exctype
    value
    traceback

    Returns
    -------

    """

    logging.exception(''.join(format_exception(exctype, value, traceback)))
    logging.warn(('#'*30))

    if isinstance(exctype, type(LamaDataException)):
        logging.warn('Lama encountered a problem with reading or interpreting some data. Please check the log files')
    else:
        logging.warn('Lama encountered an unknown problem. Please check the log files')
    sys.exit()


def command_line_agrs():
    return ' '.join(sys.argv)


def touch(file_: Path):
    if not file_.is_file():
        open(file_, 'a').close()


class LoadImage(object):
    """
    Wrapper around sitk.ReadImage which does some error checking. Takes a str or a Path
    """
    def __init__(self, img_path: Union[str, Path]):
        self.img_path = str(img_path)
        self.error_msg = None
        self.img = None
        self._read()

    def __bool__(self):
        """
        Overload this so we can do simple 'is LoadImage' to check if img loaded
        """
        if self.img is None:
            return False
        else:
            return True

    @property
    def array(self) -> np.ndarray:
        return sitk.GetArrayFromImage(self.img)

    @property
    def itkimg(self) -> sitk.Image:
        return self.img

    @property
    def direction(self) -> Tuple:
        return self.img.GetDirection()

    def _read(self):

        if self.img_path.endswith('.mnc'):
            array = read_minc.mincstats_to_numpy(self.img_path)
            self.img = sitk.GetImageFromArray(array) # temp

        if os.path.isfile(self.img_path):
            try:
                self.img = sitk.ReadImage(self.img_path)
            except RuntimeError:
                self.error_msg = "possibly corrupted file {}".format(self.img_path)

        else:
            self.error_msg = "path does not exist: {}".format(self.img_path)
            raise FileNotFoundError(f'cannot read {self.img_path}')


def check_labels_file(labelmap_file, label_info_file):
    """
    Given a label map and label file, check that the labels are all present in the info file and vice versa

    Parameters
    ----------
    labelmap_file: str
        path to labelmap
    label_info_file: str
    path to label info file

    Returns
    -------
    bool:
        True if they match
        False if not

    """

    label_map = sitk.GetArrayFromImage(sitk.ReadImage(labelmap_file))
    label_df = pd.read_csv(label_info_file, index_col=0)

    info_labels = [int(x) for x in label_df[['label']].values]
    map_labels = np.unique(label_map)

    if set(info_labels) == set(map_labels):
        return True
    else:
        return False


def plot_swarm():
    pass


def pad(array, dims):
    """
    TODO: Finish
    Pad a 3D array to dims
    Parameters
    ----------
    array: numpy.ndarray
        Array to pad
    dims: tuple (z, y, x)

    Returns
    -------
    numpy.ndarray
        padded array

    """
    p = [divmod(t, s) for s, t in zip(array.shape, dims)]

    padded = np.pad(array, ())


def write_array(array: np.ndarray, path: Union[str, Path], compressed=True, ras=True):
    """
    Write a numpy array to and image file using SimpleITK.
    If an RAS nrrd has been read by sitk, converted to numpy then read back as sitk Image it will be written out
    with the incorrect header, so it will need to be correct the directions to account for it
    """
    path = str(path)
    img = sitk.GetImageFromArray(array)
    if ras:
        img.SetDirection((-1, 0, 0, 0, -1, 0, 0, 0, 1))
    sitk.WriteImage(img, path, compressed)


def read_array( path: Union[str, Path]):
    """
    """
    path = str(path)
    return sitk.GetArrayFromImage(sitk.ReadImage(path))


def img_path_to_array(img_path: Union[str, Path]):
    if os.path.isfile(img_path):
        try:
            img = sitk.ReadImage(str(img_path))

        except RuntimeError as e:
            raise OSError(e)
        else:
            array = sitk.GetArrayFromImage(img)
            return array
    else:
        return None


def git_log() -> str:
    """
    Get the git branch and commit for logging.
    This info is in the current_commit file, which is written to after each commit using opst-commit hook.

    Returns
    -------
    the git branch, commit, and message
    """
    this_dir = Path(__file__).parent.resolve()
    git_msg_file = this_dir / 'current_commit'

    try:
        msg = ''
        with open(git_msg_file, 'r') as fh:
            for line in fh:
                msg += line
    except OSError:
        msg = f'Cannot determine git commit'

    return msg



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

    logzero.logfile(logpath)


def format_timedelta(time_delta):
    """
    Convert a datetime.timedelta to str format
    example output: "0h:0m:18s"

    Parameters
    ----------
    time_delta: datatime.timedelta

    Returns
    -------
    str: formatted timedelta

    """
    d = {}
    d["hours"], rem = divmod(time_delta.seconds, 3600)
    d["minutes"], d["seconds"] = divmod(rem, 60)
    return "{hours}h:{minutes}m:{seconds}s".format(**d)


def load_label_map_names(organ_names_path, include_terms=False):
    """
    Lod the label names csv. Remove 0 (clear label if present)
    example file format

        label_name,term
        1,l1,emap:1
        2,l2,emap:2
        3,l3,emap:4
        4,l4,emap:5
        5,l5,emap:6


    Returns
    -------
    pandas data frame
        label_num,label_name,emapa_term(optional)

    TODO: what to do if no ontology terms are available?

    """
    import pandas as pd
    df = pd.read_csv(organ_names_path)

    # Drop clear label, if present
    if df.iloc[0].label == 0:
        df.drop(0, inplace=True)

    # Check required columns are present
    required_columns = ['label', 'label_name', 'term']

    if not all(x in df for x in required_columns):
        raise ValueError(
            "The following columns are required in the label_names csv\n{}".format('\n'.join(required_columns))
        )
    return df


def mkdir_force(dir_: Union[str, Path]):
    dir_ = Path(dir_)
    shutil.rmtree(dir_, ignore_errors=True)
    dir_.mkdir(parents=True)


def mkdir_if_not_exists(dir_: Union[str, Path]):
    dir_ = Path(dir_)
    if not Path(dir_).is_dir():
        dir_.mkdir(parents=True)

def get_images_ignore_elx_itermediates(*args, **kwargs) -> Union[List[str], List[Path]]:
    return get_file_paths(*args, ignore_folders=[RESOLUTION_IMGS_DIR, IMG_PYRAMID_DIR], **kwargs)


def get_file_paths(folder: Union[str, Path], extension_tuple=('.nrrd', '.tiff', '.tif', '.nii', '.bmp', 'jpg', 'mnc', 'vtk', 'bin', 'npy'),
                   pattern: str = None, ignore_folders: Union[List, str] = []) -> Union[List[str], List[Path]]:
    """
    Given a directory return all image paths within all sibdirectories.

    Parameters
    ----------
    folder
        Where to look for images
    extension_tuple
        Select only images with these extensions
    pattern
        Do a simple `pattern in filename` filter on filenames
    ignore_folders
        do not look in folder with these names

    Notes
    -----
    Lama is currently using a mixture of Paths or str to represent filepaths. Will move all to Path.
    For now, return same type as folder input

    Do not include hidden filenames
    """
    paths = []

    if isinstance(ignore_folders, str):
        ignore_folders = [ignore_folders]


    for root, subfolders, files in os.walk(folder):

        for f in ignore_folders:
            if f in subfolders:
                subfolders.remove(f)

        for filename in files:

            if filename.lower().endswith(extension_tuple) and not filename.startswith('.'):

                if pattern:

                    if pattern and pattern not in filename:
                        continue

                paths.append(os.path.abspath(os.path.join(root, filename)))

    if isinstance(folder, str):
        return paths
    else:
        return [Path(x) for x in paths]


def check_config_entry_path(dict_, key):
    try:
        value = dict_[key]
    except KeyError:
        raise Exception("'{} not specified in config file".format(key))
    else:
        if not os.path.isdir(value):
            raise OSError("{} is not a correct directory".format(value))


def getfile_startswith(dir_: Path, prefix: str) -> Path:
    """
    Get file from a folder with a given prefix.

    Parameters
    ----------
    dir_
        Folder to search
    prefix
        The prefix to match

    Returns
    -------
    Th found file path

    Raises
    ------
    TODO: Raise exception if more than one match exists

    """
    try:
        return [x for x in dir_.iterdir() if x.name.startswith(prefix)][0]
    except IndexError as e:
        raise FileNotFoundError(f'cannot find path file starting with {prefix} in {dir_}') from e


def getfile_endswith(dir_: Path, suffix: str):
    try:
        return [x for x in dir_.iterdir() if x.name.endswith(suffix)][0]
    except IndexError as e:
        raise FileNotFoundError(f'cannot find path file ending with {suffix} in {dir_}') from e


def getfile_startswith_endswith(dir_: Path, prefix: str, suffix: str):
    try:
        return [x for x in dir_.iterdir() if x.name.endswith(suffix) and x.name.startswith(prefix)][0]
    except IndexError as e:
        raise FileNotFoundError(f'cannot find path starting with {prefix} and ending with {suffix} in {dir_}') from e



def get_inputs_from_file_list(file_list_path, config_dir):
    """
    Gte the input files

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
                raise LamaDataException

            base = line.strip()
            root_path_dict[root].append(base)
            i = 0
    for root, bases in list(root_path_dict.items()):

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


def average(img_paths: List[Path]) -> sitk.Image:
    """
    Make a mean itensity volume given a list of volume paths

    Returns
    -------
    Mean volume

    """
    img_paths = list(map(str, img_paths))

    # sum all images together
    first = LoadImage(img_paths[0])

    # Get the direction from the first image.
    direction_cos = first.direction

    summed = sitk.GetArrayFromImage(sitk.ReadImage(img_paths[0]))

    for image in img_paths[1:]:  # Ommit the first as that is in 'summed'
        np_array = sitk.GetArrayFromImage(sitk.ReadImage(image))
        try:
            summed += np_array
        except ValueError as e:
            print(("Numpy can't average this volume {0}".format(image)))

    # Now make average
    summed //= len(img_paths)
    avg_img = sitk.GetImageFromArray(summed)
    avg_img.SetDirection(direction_cos)

    return avg_img

#
# def rebuid_subsamlped_output(array, shape, chunk_size):
#     """
#
#     Parameters
#     ----------
#     array: numpy.ndarray
#         the subsampled array to rebuild
#     shape: tuple
#         the shape of the final result
#     chunk_size: int
#         the original subsampling factor
#
#     Returns
#     -------un
#     np.ndarray
#         rebuilt array of the same size of the original inputs data
#
#     """
#     out_array = np.zeros(shape)
#     i = 0
#     for z in range(0, shape[0] - chunk_size, chunk_size):
#         for y in range(0, shape[1] - chunk_size, chunk_size):
#             for x in range(0, shape[2] - chunk_size, chunk_size):
#                 out_array[z: z + chunk_size, y: y + chunk_size, x: x + chunk_size] = array[i]
#                 i += 1
#
#     return out_array

def rebuild_subsamlped_output(subsampled_array, output_array, chunk_size, mask):
    """

    Parameters
    ----------
    output_array: np.ndarray
        The 3D output array. Modified inplace
    chunk_size: size of the chunks to rebuild
    mask: np.ndarray

    """

    shape = output_array.shape
    for i, slice_ in enumerate(iterate_chunks(shape, chunk_size)):
        if not np.any(mask[slice_]):  # only mask
            output_array[slice_] = 0
        else:
            output_array[slice_] = subsampled_array[i]


def get_chunks(array, chunk_size, mask):
    """
    Get an iterator of chunks of data from array

    Parameters
    ----------
    array: np.ndarray
        array to chunck or to rebuild
    chunk_size: int
    mask: np.ndarray

    Returns
    -------
    iterator<ndarray>
    """
    shape = array.shape

    for slice_ in iterate_chunks(shape, chunk_size):

        if not np.any(mask[slice_]):
            continue

        else:
            yield array[slice_]


def iterate_chunks(shape, chunk_size):
    for z in range(0, shape[0] - chunk_size, chunk_size):
        for y in range(0, shape[1] - chunk_size, chunk_size):
            for x in range(0, shape[2] - chunk_size, chunk_size):
                slice_ = np.s_[z: z + chunk_size, y: y + chunk_size, x: x + chunk_size]
                yield slice_


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
        for root, basenames in root_names_dict.items():
            fh.write({'dir:{}\n'.format(root)})
            for base in basenames:
                fh.write('{}\n'.format(base))


def csv_read_lines(path):
    """
    Read lines from a csv. Each line is assumed to have one entry only, such as a specimen id
    ----------
    path: str
        path to csv file

    Returns
    -------
    list of lines
    """
    lines = []
    with open(path, 'r') as fh:
        for line in fh:
            lines.append(line.strip())
    return lines


def csv_to_pandas(path):
    """
    Pandas cannot open csv files that have locks on them. For example if they are currently opened by LibreOffice
    This is a workaround

    Parameters
    ----------
    path: str
        path to csv

    Returns
    -------
    pandas.Dataframe

    Not yet finished

    """
    import pandas as pd
    with open(path, 'r') as fh:
        try:
            df = pd.read_csv(fh)
        except:
            pass


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
    with open(path, 'r') as fh:
        reader = csv.reader(fh)
        for line in reader:
            if not line:
                continue
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
        stripped = file_name.rstrip(ext)
        return stripped
    else:
        return file_name


def test_installation(app):
    try:
        sub.check_output([app])
    except Exception as e:  # can't seem to log CalledProcessError
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


def service_shutdown(signum, frame):
    """
    Catches termination signal, writes to log.

    Parameters
    ----------
    signum
    frame

    Returns
    -------

    """
    logging.info(f"Caught {signal.Signals(signum)} signal\n. Lama exiting")
    raise SystemExit


def available_memory() -> float:
    """
    Get the number of memory availbe in bytes
    """
    return psutil.virtual_memory().available


def download_and_extract_zip(url: Path, out_dir: Path):

    print('Downloading data')
    
    
    http.HTTPConnection._http_vsn = 10
    http.HTTPConnection._http_vsn_str = 'HTTP/1.0'
    
    remotezip = urllib.request.urlopen(url)
    
    zipinmemory = io.BytesIO(remotezip.read())
    
    http.HTTPConnection._http_vsn = 11
    http.HTTPConnection._http_vsn_str = 'HTTP/1.1'    
    zip = zipfile.ZipFile(zipinmemory)
    zip.extractall(out_dir)
     
    
    print(f'Data downloaded and extracted to {out_dir}')


def is_number(value):
    """

    Parameters
    ----------
    value

    Returns
    -------

    """
    try:
        decimal.Decimal(value)
    except decimal.DecimalException:
        return False
    else:
        return True


def truncate_str(string: str, max_len: int) -> str:
    string = str(string)
    if len(string) > max_len:
        string = f'{string[:max_len - 4]} .. {string[len(string) - 4: ]}'
    return string


def cfg_load(cfg) -> Dict:
    """
    There are 2 types of config file used in the project yaml an toml. Will move to al tml at some point

    This function wraps around both

    Returns
    -------
    Dictionary config
    """
    cfg = Path(cfg)

    if not cfg.is_file():
        raise FileNotFoundError(f'Cannot find required config file: {cfg}')

    if Path(cfg).suffix == '.yaml':

        try:
            with open(cfg, 'r') as fh:
                return yaml.load(fh, Loader=yaml.FullLoader)
        except Exception as e:
            raise ValueError("can't read the config file - {}".format(e))

    elif Path(cfg).suffix == '.toml':
        try:
            return toml.load(cfg)
        except Exception as e:
            raise ValueError("can't read the config file - {}".format(e))

    else:
        raise ValueError('Config file should end in .toml or .yaml')


def bytesToGb(numberOfBytes, precision = 3):
    return round(numberOfBytes / (1024 ** 3), precision)
    
def logMemoryUsageInfo():
    proc = psutil.Process(os.getpid())
    
    procMemInfo = proc.memory_info()
    globalMemInfo = psutil.virtual_memory()

    totalMem = bytesToGb(globalMemInfo.total, 5)
    totalMemAvailable = bytesToGb(globalMemInfo.available, 5)
    procMemRSS = bytesToGb(procMemInfo.rss, 5)
    procMemVMS = bytesToGb(procMemInfo.vms, 5)
    procMemData = bytesToGb(procMemInfo.data, 5)
    
    logging.info(f"Memory: Process Resident: {procMemRSS}, Process Virtual: {procMemVMS}, Process data: {procMemData}.  Total Available: {totalMemAvailable}")
