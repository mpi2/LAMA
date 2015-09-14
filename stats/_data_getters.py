
import sys
import os
import SimpleITK as sitk
import numpy as np
# Hack. Relative package imports won't work if this module is run as __main__
sys.path.insert(0, os.path.abspath('..'))
import common


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
        self.wt_data_dir = wt_data_dir
        self.mut_data_dir = mut_data_dir
        self.wt_paths, self.mut_paths = self.get_data_paths()
        self.wt_data, self.mut_data = self.generate_data()

    def get_data_paths(self):
        """
        Get paths to the data
        """
        wt_paths = common.GetFilePaths(self.wt_data_dir)
        mut_paths = common.GetFilePaths(self.mut_data_dir)
        return wt_paths, mut_paths

    def generate_data(self):
        """
        Process the data so it's in scalar form to be analysed by the stats classes
        """
        wt_data = []
        mut_data = []
        for wt_data_path in self.wt_paths():
            wt_data.append(self.get_data(wt_data_path))
        for mut_data_path in self.mut_paths:
            mut_data.append(self.get_data(mut_data_path))

        return wt_data, mut_data

    def get_data(self, path_):
        raise NotImplementedError

    @staticmethod
    def blur_volume(img):
        """
        Create memory-mapped volumes
        """
        # previous: 1.0, 8, 0.001
        blurred = sitk.DiscreteGaussian(img, 0.5, 4, 0.01, False)
        return blurred


class ScalarDataGetter(AbstractDataGetter):
    """
    Processess the image intensity data:
        - The normnalised registered images
        - Spatial jacobians
    """
    def __init__(self, *args):
        super(ScalarDataGetter, self).__init__(*args)

    def get_data(self, path_):
        """
        For Intensity and jacobians, we already have scalr data that can go into the stats calculations as-is
        So just return it
        """
        return self.blur(sitk.GetArrayFromImage(sitk.ReadImage(path_)))


class DeformationDataGetter(AbstractDataGetter):
    """
    Process the deformations fields generated during registration
    """
    def __init__(self, *args):
        super(DeformationDataGetter, self).__init__(*args)

    def get_data(self, path_):
        """
        Calculates the deformation vector magnitude at each voxel position
        """
        arr = sitk.GetArrayFromImage(sitk.ReadImage(path_))
        return np.sqrt((arr*arr).sum(axis=3))


class GlcmDataGetter(AbstractDataGetter):
    """
    Get data from the grey level co-occurence matrices
    """
    def __init__(self, *args):
        super(GlcmDataGetter, self).__init__(*args)

    def get_data(self, path_):
        """
        The data will be in this form:
        A numpy npz for both wt and mutant. Aeach containing glcms for each block for each specimen
        """
        pass # Need to call the glcm module and extract the features



