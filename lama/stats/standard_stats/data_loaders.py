"""
Load data for the respective datatypes

Currently: Converts output from registration to 8bit. As the registration pipeline only accepts 8bit images at the moment
this is ok. When we change to allow 16 bit images, we may have to change a few things in here
"""

import os
from pathlib import Path
from typing import Union, List, Iterator, Tuple

import numpy as np
from addict import Dict
import SimpleITK as sitk
from logzero import logger as logging
import pandas as pd
import yaml

from lama import common
from lama.img_processing.normalise import normalise
from lama.img_processing.misc import blur
from lama.paths import iterate_over_specimens


GLCM_FILE_SUFFIX = '.npz'
DEFAULT_FWHM = 100  # um
DEFAULT_VOXEL_SIZE = 14.0
IGNORE_FOLDER = 'resolution_images'

# This is temporary. need ot centralise all path stuff in lama.paths


class InputData:
    """
    Holds the input data that will be analysed.
    Just a wrpper around a pandas DataFrame with methods to get various elements
    """
    def __init__(self, data: np.ndarray,
                 info: pd.DataFrame,
                 line: str,
                 shape: Tuple):
        """
        Holds the input data to be used in the stats tests
        Parameters
        ----------
        data
            2D array
                row: specimens
                columns: data points
        info
            columns:
                - specimen (index)
                - staging (the staging metric)
                - line
                - genotype


        """
        self.data = data
        self.info = info
        self.shape = shape
        self.line = line

        if data.shape[0] != len(info):
            raise ValueError

    def mutant_ids(self):
        return self.info[self.info.genotype == 'mutant'].index

    def genotypes(self):
        return self.info.genotype

    def chunks(self, chunk_size) -> Iterator[np.ndarray]:
        """
        Split the return the data in chunks

        Yields
        -------
        Chunks split column-wise (axis=1)

        Notes
        -----
        np.array_split return a view on the data not a copy
        """
        num_pixels = self.data.shape[1]
        num_chunks = num_pixels / chunk_size if num_pixels > chunk_size else 1
        chunks = np.array_split(self.data, num_chunks, axis=1)

        for data_chunk in chunks:
            yield data_chunk


class DataLoader:
    """
    Parent class for loading in data

    Notes
    -----
    TODO: add support for memory mapping data
    """
    def __init__(self,
                 wt_dir: Path,
                 mut_dir: Path,
                 mask: np.ndarray,
                 config: Dict):

        self.wt_dir = wt_dir
        self.mut_dir = mut_dir
        self.config = config
        self.mask = mask  # 3D mask
        self.shape = None

        self.blur_fwhm = config.get('blur', DEFAULT_FWHM)
        self.voxel_size = config.get('voxel_size', DEFAULT_VOXEL_SIZE)

    @staticmethod
    def factory(type_):
        if type_ == 'intensity':
            return IntensityDataGetter
        elif type_ == 'jacobians':
            return JacobianDataGetter
        elif type_ == 'organ_volume':
            raise NotImplementedError('Organ volumes data getter not working yet')

    def line_iterator(self) -> Iterator[InputData]:
        """

        Returns
        -------

        """
        raise NotImplementedError

    def _read(self, paths) -> np.ndarray:
        """
        - Read in the voxel-based data into 3D arrays
        - Apply guassian blur to the 3D image
        - mask
        - Unravel


        Parameters
        ----------
        paths
            Path to load

        Returns
        -------
        2D array
            columns: Data points
            rows: specimens
        """

        images = []

        for data_path in paths:
            loader = common.LoadImage(data_path)

            if not self.shape:
                self.shape = loader.array.shape

            blurred_array = blur(loader.array, self.blur_fwhm, self.voxel_size)
            masked = blurred_array[self.mask != False]

            images.append(masked.ravel())

        return np.array(images)


class JacobianDataGetter(DataLoader):
    """
    Process the Spatial Jacobians generated during registration
    """
    def __init__(self, *args):
        super(JacobianDataGetter, self).__init__(*args)
        self.datatype = 'jacobians'

    def line_iterator(self) -> InputData:

        # Get the wild type paths and ids
        wt_info = get_paths(self.wt_dir, self.datatype, str(self.config['jac_folder']))
        masked_wt_data = self._read(wt_info['path'])

        # Get all the mutant paths and ids
        mut_info = get_paths(self.mut_dir, self.datatype, str(self.config['jac_folder']))

        # Iterate over the lines
        mut_gb = mut_info.groupby('line')
        for line, mut_df in mut_gb:
            masked_mut_data = self._read(mut_df['path'])

            # Make dataframe of specimen_id, genotype, staging
            wt_staging = get_staging_data(self.wt_dir)
            wt_staging['genotype'] = 'wildtype'
            mut_staging = get_staging_data(self.mut_dir, line=line)
            mut_staging['genotype'] = 'mutant'

            info = pd.concat((wt_staging,  mut_staging))
            # Id there is a value column, change to staging. TODO: make lama spitout staging header instead of value
            if 'value' in info:
                info.rename(columns={'value': 'staging'}, inplace=True)

            data = np.vstack((masked_wt_data, masked_mut_data))
            input_ = InputData(data, info, line, self.shape)
            yield input_


















class IntensityDataGetter(DataLoader):
    """
    Processess the image intensity data:
        - The normnalised registered images
    """
    def __init__(self, *args):
        super(IntensityDataGetter, self).__init__(*args)

    def _get_normalised_data(self):
        """
        Normalise both WT and mut data to the same roi and memory map the results
        Returns
        -------
        dirs: tuple
            [0] normalised wt dir
            [1] normlaised mut dir
        """
        wt_norm_paths, mut_norm_paths = normalise(self.wt_paths, self.mut_paths,
                                                  self.normalisation_dir, self.normalisation_roi)
        return wt_norm_paths, mut_norm_paths

    def load_data(self, wt_dir: Path, mut_dir: Path) -> InputData:
        """
        Load the data into an InputData data object and return

        For Intensity and jacobians, we already have scalar data
        that can go into the stats calculations as is with just some blurring

        Returns
        -------
        The input data values and associated groups information
        """

        # loop over wt and then mutant paths


        if self.normalisation_roi is not None:
            wt_paths, mut_paths = self._get_normalised_data()

        logging.info("\nBlurring and masking data\nVoxel size (um):{}\nblur radius fwhm (um): {}\n".format(
            self.voxel_size, self.blur_fwhm))
        masked_wt_data = load(wt_paths)
        masked_mut_data = load(mut_paths)

        loader = common.LoadImage(wt_paths[0])
        if not loader:
            logging.error("Problem getting data for stats: {}".format(loader.error_msg))
        self.shape = loader.array.shape

        return masked_wt_data, masked_mut_data




class DeformationDataGetter(DataLoader):
    """
    Process the deformations fields generated during registration
    """
    def __init__(self, *args):
        super(DeformationDataGetter, self).__init__(*args)

    def _get_data(self, wt_paths, mut_paths):
        """
        Calculates the deformation vector magnitude at each voxel position
        """
        def load(paths):
            array = []
            for data_path in paths:
                arr_32bit = common.img_path_to_array(data_path).astype(np.float32)
                if arr_32bit.ndim != 4:
                    msg = "The deformation files are not 4D. Is the stats config correct?"
                    logging.error(msg)
                    raise ValueError(msg)
                vector_magnitudes = np.sqrt((arr_32bit*arr_32bit).sum(axis=3))
                v_img = sitk.GetImageFromArray(vector_magnitudes)
                blurred_array = self._blur_volume(v_img).ravel()
                masked = blurred_array[self.mask != False]
                memmap_array = self._memmap_array(masked)
                array.append(memmap_array)
            return array

        masked_wt_data = load(wt_paths)
        masked_mut_data = load(mut_paths)
        self.shape = common.img_path_to_array(wt_paths[0]).shape[0:3]  # 4th dimension is the deformation vector
        return masked_wt_data, masked_mut_data


class GlcmDataGetter(DataLoader):
    """
    Get data from the grey level co-occurence matrices
    """
    def __init__(self, *args):
        super(GlcmDataGetter, self).__init__(*args)

    def _get_data(self, wt_paths, mut_paths):
        # Get the glcm config from the data directories
        dir_ = os.path.split(wt_paths[0])[0]
        glcm_metadata_path = os.path.join(dir_, 'glcm.yaml')
        with open(glcm_metadata_path, 'r') as fh:
            glcm_metadata = yaml.load(fh)
        shape = glcm_metadata['original_shape']
        self.shape = shape

        def load(paths):

            result = []

            for data_path in paths:

                glcm_features = np.fromfile(data_path, dtype=np.float32)
                result.append(glcm_features.ravel())

            return result

        wt_data = load(wt_paths)
        mut_data = load(mut_paths)

        return wt_data, mut_data


def get_paths(root_dir: Path, data_type: str, subfolder: str) -> pd.DataFrame:
    """
    Get the data paths for the data type specified by 'datatype'

    Parameters
    ----------
    root_dir
        Registration output directory to search
    datatype
        what data to look for

    Returns
    -------
    The Paths to the data

    Raises
    ------
    FileNotFoundError if any data is missing
    """
    reg_out_dir = root_dir / 'output'
    specimen_info = []

    for line_dir in reg_out_dir.iterdir():

        if not line_dir.is_dir():
            continue

        for spec_dir in line_dir.iterdir():

            if str(spec_dir).endswith('_'):  # previous 'stats_' directory
                continue

            spec_out_dir = spec_dir / 'output'

            if not spec_out_dir.is_dir():
                raise FileNotFoundError(f'Cannot find output directory for {spec_dir}')

            data_dir = spec_out_dir / data_type / subfolder
            if not data_dir.is_dir():
                raise FileNotFoundError(f'Cannot find data directory: {data_dir}')

            # Data file will have same name as specimen
            data_file = data_dir.glob(f'{spec_dir.name}*')

            if data_file:
                # For each specimen we have: id, line and the data file path
                specimen_info.append([spec_dir.name, line_dir.name, data_file.__next__()])

            else:
                raise FileNotFoundError(f'Data file missing: {data_file}')

    df = pd.DataFrame.from_records(specimen_info, columns=['specimen', 'line', 'path'])
    return df


def load_mask(parent_dir, mask_path: Path) -> np.ndarray:
    """
    Mask is used in multiple datagetter so weload it independently of the classes.

    Raises
    ------
    ValueError if mask contains anything other than ones and zeroes

    Returns
    -------
    mask 3D
    """
    mask = common.LoadImage(parent_dir / mask_path).array

    if set([0, 1]) != set(np.unique(mask)):
        logging.error("Mask image should contain only ones and zeros ")
        raise ValueError("Mask image should contain only ones and zeros ")

    return mask


def get_staging_data(root: Path, line=None) -> pd.DataFrame:
    """
    Collate all the staging data from a folder.

    root
        The root directory to search
    line
        Only select staging data for this line

    """

    output_dir = root / 'output'

    dataframes = []

    for line_dir, specimen_dir in iterate_over_specimens(output_dir):
        if line and line_dir.name != line:
            continue

        staging_file = specimen_dir / 'output' / common.STAGING_INFO_FILENAME

        if not staging_file.is_file():
            raise FileNotFoundError(f'Cannot find organ volume file {staging_file}')

        df = pd.read_csv(staging_file, index_col=0)
        df['line'] = line_dir.name
        dataframes.append(df)

    # Write the concatenated organ vol file to single csv
    staging = pd.concat(dataframes)

    outpath = output_dir / common.ORGAN_VOLUME_CSV_FILE
    staging.to_csv(outpath)

    return staging


def Gamma2sigma(Gamma):
    '''Function to convert FWHM (Gamma) to standard deviation (sigma)'''
    return Gamma * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)