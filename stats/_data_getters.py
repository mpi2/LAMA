
import sys
import os
import SimpleITK as sitk
import numpy as np
# Hack. Relative package imports won't work if this module is run as __main__
sys.path.insert(0, os.path.abspath('..'))
import common
from glcm3d import ContrastTexture
from os.path import join


GLCM_FILE_SUFFIX = '.npz'


class AbstractDataGetter(object):
    """
    Gets the data. Could be scalar, vector or texture
    """
    def __init__(self, wt_data_dir, mut_data_dir):
        """
        Parameters
        ----------
        wt_mut_data_dir: str
            path to folder containing data volumes
        mut_data_dir: str
            path to folder containing data volumes

        Notes
        _____
        Data can exist in sub-folders
        """
        self.shape = None
        self.wt_data_dir = wt_data_dir
        self.mut_data_dir = mut_data_dir
        self.wt_paths, self.mut_paths = self._get_data_paths()
        self.wt_data, self.mut_data = self._generate_data()


    def _get_data_paths(self):
        """
        Get paths to the data
        """
        wt_paths = common.GetFilePaths(self.wt_data_dir)
        mut_paths = common.GetFilePaths(self.mut_data_dir)
        return wt_paths, mut_paths

    def _generate_data(self):
        """
        Process the data so it's in scalar form to be analysed by the stats classes

        Returns
        -------
        mut and wt data are in lists. each specimen data file should be 3d reshaped
        """
        wt_data = self._get_data(self.wt_paths)
        mut_data = self._get_data(self.mut_paths)


        # for wt_data_path in self.wt_paths:
        #     wt_data.append()
        # for mut_data_path in self.mut_paths:
        #     mut_data.append(self._get_data(mut_data_path))

        return wt_data, mut_data

    @property
    def wildtype_data(self):
        return self.wt_data

    @property
    def mutant_data(self):
        return self.mut_data

    def _get_data(self, path_):
        """
        This is the main method to overide. Convert the path to some form of data
        """
        raise NotImplementedError

    @staticmethod
    def _blur_volume(img):
        """

        """
        # previous: 1.0, 8, 0.001
        blurred = sitk.DiscreteGaussian(img, 0.5, 4, 0.01, False)
        return sitk.GetArrayFromImage(blurred)


class ScalarDataGetter(AbstractDataGetter):
    """
    Processess the image intensity data:
        - The normnalised registered images
        - Spatial jacobians
    """
    def __init__(self, *args):
        super(ScalarDataGetter, self).__init__(*args)

    def _get_data(self, paths):
        """
        For Intensity and jacobians, we already have scalr data that can go into the stats calculations as-is
        So just return it

        Returns
        -------
        blurred_array: np ndarry
            bluured image array
        """
        result = []
        self.shape = common.img_path_to_array(paths[0]).shape
        for data_path in paths:
            blurred_array = self._blur_volume(sitk.ReadImage(data_path))
            result.append(blurred_array)
        return result


class DeformationDataGetter(AbstractDataGetter):
    """
    Process the deformations fields generated during registration
    """
    def __init__(self, *args):
        super(DeformationDataGetter, self).__init__(*args)

    def _get_data(self, paths):
        """
        Calculates the deformation vector magnitude at each voxel position
        """
        self.shape = common.img_path_to_array(paths[0]).shape[0:3]  # 4th dimension is the deformation vector
        result = []
        for data_path in paths:
            arr = common.img_path_to_array(data_path)
            vector_magnitudes = np.sqrt((arr*arr).sum(axis=3))
            blurred_array = self._blur_volume(sitk.GetImageFromArray(vector_magnitudes))
            result.append(blurred_array)
        return result


class GlcmDataGetter(AbstractDataGetter):
    """
    Get data from the grey level co-occurence matrices
    """
    def __init__(self, *args):
        super(GlcmDataGetter, self).__init__(*args)

    def _get_data(self, paths):
        """
        Parameters
        ----------
        path_: str
            this will be a path to a GLCM array containing multiple specimen GLCMs
        """
        # For now, just use the contrast measure
        first = True
        npz = np.load(paths[0]) # For glcm analysis, there's only one data file in the list
        header = npz['header'][()]
        self.glcm_chunk_size = header['chunk_size']
        contrast_measures = []
        if first:
            self.shape = header['image_shape']
            first = False
        if self.shape != header['image_shape']:
            raise ValueError('"Shapes of volumes are not the same. They need to be the same for stats analysis"')
        for specimen in npz['data']:
            contrast_measures.append(ContrastTexture(specimen, header).get_results())
        return np.array(contrast_measures)

    def _get_data_paths(self):
        """
        Override this as we don't want to get a list of file, we jsut want the .npz glcm array
        """
        wt_glcms = [join(self.wt_data_dir, x) for x in os.listdir(self.wt_data_dir) if x.endswith(GLCM_FILE_SUFFIX)]
        if len(wt_glcms) != 1:
            print "There should only be one GLCM array in the folder"
            return
        mut_glcms = [join(self.mut_data_dir, x) for x in os.listdir(self.mut_data_dir) if x.endswith(GLCM_FILE_SUFFIX)]
        if len(mut_glcms) != 1:
            print "There should only be one GLCM array in the folder"
            return
        else:
            return wt_glcms, mut_glcms




