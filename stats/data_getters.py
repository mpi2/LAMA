"""
Currently: Converts output from registration to 8bit. As the registration pipeline only accepts 8bit images at the moment
this is ok. When we change to allow 16 bit images, we may have to change a few things in here
"""
import logging
import os
import sys
import tempfile
import datetime
import SimpleITK as sitk
import numpy as np
import scipy.ndimage as ndimage

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

# Hack. Relative package imports won't work if this module is run as __main__
sys.path.insert(0, os.path.abspath('..'))
import common
from img_processing.normalise import normalise

GLCM_FILE_SUFFIX = '.npz'
DEFAULT_FWHM = 100  # 100 um
DEFAULT_VOXEL_SIZE = 28.0
IGNORE_FOLDER = 'resolution_images'



class AbstractDataGetter(object):
    """
    Gets the data. Could be scalar, vector or texture
    """
    def __init__(self, wt_data_paths, mut_data_paths, mask, volorder=None, voxel_size=None,
                 subsampled_mask=None, subsample_int=None, blur_fwhm=None, memmap_dir=None):
        """
        Parameters
        ----------
        wt_data_paths: list
            paths to the wild type data files
        mut_data_paths: list
            paths to the wild type data files
        mask: np.ndarray
            flattened mask
        volorder: list/None  -> do we need this?
            if list, reorder the data once got, so it's the same order as the groups file (for linear models etc)
        subsample: bool/int
             False: do not subsample data
             int: subsampling size. subsample and additionaly provide this along with unsubsample data
        memmap_dir: str
            Path to folder to create temporary file for numpy array memmapping. Uses /temp if not set
        """
        if blur_fwhm:
            self.blur_fwhm = blur_fwhm
        else:
            self.blur_fwhm = DEFAULT_FWHM
        if voxel_size:
            self.voxel_size = voxel_size
        else:
            self.voxel_size = DEFAULT_VOXEL_SIZE
        self.mask = mask
        if self.mask.max() != 1:
            logging.error("Mask image should contain only ones and zeros ")
            raise ValueError("Mask image should contain only ones and zeros ")
        # if self.mask.min() != 0:
        #     logging.error("Mask image should contain only ones and zeros ")
        #     sys.exit()

        self.subsample_int = subsample_int
        self.subsampled_mask = subsampled_mask
        self.volorder = volorder
        self.shape = None
        self.subsampled_shape = None
        self.wt_paths = wt_data_paths
        self.mut_paths = mut_data_paths
        self.voxel_size = voxel_size
        self.memmap_dir = memmap_dir

    def set_data(self):
        """
        Read in all the data into numpy arrays.
        Order the data such that it's the same as the groups file that is given to R
        Returns
        -------

        """
        self._generate_data()

        # Check if numpy of paths == volumes listed in groups.csv. If volorder == None, we don't have a groups.csv file

        if self.volorder:
            total_volumes = len(self.masked_mut_data) + len(self.masked_wt_data)

            if total_volumes != len(self.volorder):
                logging.error("Number of data files found is not the same as specified in groups file\nAre the names\
                                  in the file correct!\nnumber of volumes:{}  number of files in csv{}".format(
                    total_volumes, len(self.volorder)))

                logging.info("wt vols: {}\nmut vols: {}\ngroups file entries".format(
                    "\n".join(self.wt_paths), "\n".join(self.mut_paths), "\n".join(self.volorder)))

    def set_normalisation_roi(self, roi, normalisation_dir):
        self.normalisation_roi = roi
        self.normalisation_dir = normalisation_dir

    def reorder_paths(self, paths):
        """
        Reorder the volume paths based on the order given in groups.csv.
        If there are duplicates in groups.csv, die as we can't do the stats
        """
        ordered_paths = []
        if len(set(self.volorder)) != len(self.volorder):
            logging.error("Error! There are duplcates in the groups csv file. Exiting!")
            sys.exit()
        for vol in self.volorder:
            # Get the basename without extension as this moght change (eg from .tif to .nrrd during registration)
            v = os.path.splitext(vol)[0]
            for p in paths:
                pbase = os.path.splitext(os.path.basename(p))[0]
                # Todo: check if all files in csv are in file path list and vice versa
                if v == pbase:
                    ordered_paths.append(p)
        return ordered_paths

    def _generate_data(self):
        """
        Process the data so it's in scalar form to be analysed by the stats classes

        Returns
        -------
        mut and wt data are in lists. each specimen data file should be 3d reshaped
        """
        logging.info('reading data')
        start_time = datetime.datetime.now()
        self.masked_wt_data, self.masked_mut_data = self._get_data(self.wt_paths, self.mut_paths)
        end_time = datetime.datetime.now()
        total_seconds = (end_time - start_time).seconds
        load_mins = int(total_seconds / 60)
        if total_seconds < 1:
            load_secs = 0
        else:
            load_secs = 60 % total_seconds
        logging.info('loading of data took {} mins:{} secs'.format(load_mins, load_secs))

    def _flatten(self, arrays):
        one_d = []
        for arr in arrays:
            f = arr.flatten()  # try ravel to get a view rather than copy
            one_d.append(f)
        stacked = np.vstack(one_d)
        return stacked

    @property
    def wildtype_data(self):
        return self.wt_data

    @property
    def mutant_data(self):
        return self.mut_data


    # def _get_zscore_overlay(self):
    #     mut_mean = np.mean(self.mut_data, axis=0)
    #     wt_mean = np.mean(self.wt_data, axis=0)
    #     wt_std = np.std(self.wt_data, axis=0)
    #     zscores = (mut_mean - wt_mean) / wt_std
    #     return zscores

    def _get_data(self, path_):
        """
        This is the main method to overide. Convert the path to some form of data
        """
        raise NotImplementedError

    def _blur_volume(self, img):
        """
        https://matthew-brett.github.io/teaching/random_fields.html

        Parameters
        ----------
        img: SimpleITK Image
        """
        array = sitk.GetArrayFromImage(img)
        fwhm_in_voxels = self.blur_fwhm / self.voxel_size

        sd = fwhm_in_voxels / np.sqrt(8. * np.log(2))  # sigma for this FWHM
        blurred = ndimage.filters.gaussian_filter(array, sd, mode='constant', cval=0.0)

        return blurred

    # def _blur_volume(self, img):
    #     """
    #     Parameters
    #     ----------
    #     img: SimpleITK Image
    #     """
    #
    #     fwhm_in_voxels = self.blur_fwhm / self.voxel_size
    #
    #     sigma = Gamma2sigma(fwhm_in_voxels)
    #
    #     blurred = sitk.DiscreteGaussian(img, variance=sigma, useImageSpacing=False)
    #
    #     # return sitk.GetArrayFromImage(blurred)
    #
    #     return sitk.GetArrayFromImage(blurred)

    def _memmap_array(self, array):
        # Remove this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        return array
        t = tempfile.TemporaryFile()
        m = np.memmap(t, dtype=array.dtype, mode='w+', shape=array.shape)
        m[:] = array[:]
        return m


class IntensityDataGetter(AbstractDataGetter):
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
        print('normalising images\nWTs:\n{}\n\nMutants: {}'.format('\n'.join([x for x in self.wt_paths]),
                                                                   '\n'.join([x for x in self.mut_paths])))
        wt_norm_paths, mut_norm_paths = normalise(self.wt_paths, self.mut_paths,
                                                  self.normalisation_dir, self.normalisation_roi)
        return wt_norm_paths, mut_norm_paths

    def _get_data(self, wt_paths, mut_paths):
        """
        For Intensity and jacobians, we already have scalr data that can go into the stats calculations as-is
        So just return it

        Returns
        -------
        blurred_array: np ndarry
            bluured image array
        """
        def load(paths):
            array = []
            for data_path in paths:
                loader = common.LoadImage(data_path)
                if not loader:
                    logging.error("Problem getting data for stats: {}".format(loader.error_msg))
                data8bit = loader.img
                blurred_array = self._blur_volume(data8bit).ravel()
                masked = blurred_array[self.mask != False]
                memmap_array = self._memmap_array(masked)
                array.append(memmap_array)
            return array

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


class JacobianDataGetter(AbstractDataGetter):
    """
    Process the Spatial Jacobians generated during registration
    """
    def __init__(self, *args):
        super(JacobianDataGetter, self).__init__(*args)

    def _get_data(self,  wt_paths, mut_paths):

        def load(paths):
            array = []
            initial_vol = common.img_path_to_array(paths[0])
            if initial_vol is None:
                err_msg = "Error reading jacobian file: {}".format(paths[0])
                print err_msg
                logging.error(err_msg)
                raise common.LamaDataException()
            self.shape = initial_vol.shape
            for data_path in paths:
                img = sitk.ReadImage(data_path)
                blurred_array = self._blur_volume(img).ravel()
                masked = blurred_array[self.mask != False]
                memmap_array = self._memmap_array(masked)
                array.append(memmap_array)
            return array

        masked_wt_data = load(wt_paths)
        masked_mut_data = load(mut_paths)
        loader = common.LoadImage(wt_paths[0])
        if not loader:
            logging.error("Problem getting data for jacobian stats: {}".format(loader.error_msg))
            raise ValueError("Problem getting data for jacobian stats: {}".format(loader.error_msg))
        self.shape = loader.array.shape

        return masked_wt_data, masked_mut_data


class DeformationDataGetter(AbstractDataGetter):
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


class AngularDataGetter(AbstractDataGetter):
    """
    Process the deformations fields generated during registration
    """
    def __init__(self, *args):
        super(AngularDataGetter, self).__init__(*args)

    def _get_data(self, paths):
        """
        Calculates the deformation vector magnitude at each voxel position
        """
        result = []

        self.shape = common.img_path_to_array(paths[0]).shape[0:3]  # 4th dimension is jacobian matrix
        for data_path in paths:
            arr = common.img_path_to_array(data_path)
            a = arr
            v1 = a.take(0, axis=3).ravel()
            v2 = a.take(1, axis=3).ravel()
            angles = np.arctan2(v1, v2)
            masked = np.array(angles)[self.mask != False]
            memmap_array = self._memmap_array(masked)
            result.append(memmap_array)
        return result


class GlcmDataGetter(AbstractDataGetter):
    """
    Get data from the grey level co-occurence matrices
    """
    def __init__(self, *args):
        super(GlcmDataGetter, self).__init__(*args)

    def _get_data(self, paths):

        result = []

        for data_path in paths:
            glcm_features = np.fromfile(data_path, dtype=np.float32)
            result.append(glcm_features.ravel())
        return result


def Gamma2sigma(Gamma):
    '''Function to convert FWHM (Gamma) to standard deviation (sigma)'''
    return Gamma * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)