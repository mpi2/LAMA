#!/usr/bin/env python

import os
import sys
from collections import OrderedDict
import tempfile
from pathlib import Path
from typing import Dict, List
from itertools import accumulate

from logzero import logger as logging
import numpy as np
import SimpleITK as sitk
from lama.paths import specimen_iterator
from lama import common
from lama.registration_pipeline.validate_config import LamaConfig
from scipy import ndimage

try:
    from skimage.draw import line_aa

    skimage_available = True
except ImportError:
    skimage_available = False


class Normaliser:
    def __init__(self):

        self.mask = None

    @staticmethod
    def factory(type_, data_type: str):

        if data_type == 'intensity':

            # If passing an ROI as as list
            if isinstance(type_, (list,)):  # Not working at the moment
                if len(type_) != 3:
                    return None
                return None  # RoiNormalise

            elif type_ == 'mask':
                return IntensityMaskNormalise()
            elif type_ == 'histogram':
                return IntensityHistogramMatch()

        else:
            return None

    def memorymap_data(self, lama_root_dir: Path) -> Dict[str, np.memmap]:
        """
        Iterate over output folder getting each ...........
        Parameters
        ----------
        lama_root_dir

        Returns
        -------

        """

        imgs = OrderedDict()

        for line_dir, spec_dir in specimen_iterator(lama_root_dir):
            config_file = common.getfile_endswith('.toml')  # Get the Lama config from the specimen directory
            config = LamaConfig(config_file)
            reg_dir = config['root_reg_dir']
            basename = os.path.basename(imgpath)
            loader = common.LoadImage(imgpath)

            if not loader:
                logging.error("Problem normalising image: {}".format(loader.error_msg))
                sys.exit()
            arr = loader.array
            t = tempfile.TemporaryFile()
            m = np.memmap(t, dtype=arr.dtype, mode='w+', shape=arr.shape)
            m[:] = arr
            imgs[basename] = m
        return imgs

    def add_reference(self, ref: np.ndarray):
        """
        Add the reference data to
        Returns
        -------

        """
        raise NotImplementedError

    def normalise(self) -> np.ndarray:
        raise NotImplementedError


class NonRegMaskNormalise(Normaliser):
    """
    Normalise a set of volumes to the mean of voxel included in a mask.
    In this case each volume needs its mask
    as its not deformed.

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, *kwargs)
        self.reference_mean = None

    def add_reference(self, ref: np.ndarray, ref_mask: np.ndarray):
        """
        Add the

        Parameters
        ----------

        Returns
        -------
        """
        logging.info('normalising intensity data to mean of the mask')

        # so when we add the reference, we're not storing the image
        # so we can slice it to make computation time quicker
        img = sitk.GetArrayFromImage(ref)
        mask = sitk.GetArrayFromImage(ref_mask)

        s = ndimage.find_objects(mask)[0]

        mask = mask[s[0].start:s[0].stop,
               s[1].start:s[1].stop,
               s[2].start:s[2].stop]
        img = img[s[0].start:s[0].stop,
              s[1].start:s[1].stop,
              s[2].start:s[2].stop]

        # test if this improves speed

        # ignore vals outside of mask
        img[(mask != 1)] = 0

        self.reference_mean = np.mean(img)

    def normalise(self, volumes: List[np.ndarray], masks: List[np.ndarray],
                  fold: bool = False, temp_dir: Path = None):
        """
        given paths to registered images, apply linear normalisation so that the mean of the roi across all images are
        the same.

        Create new diretories and place the normalised images in

        Parameters
        ----------
        volumes : list of imgs
        masks: list of masks
        fold : performs fold difference if true

        Returns
        -------
        None
            Data is normalised in-place
        """

        logging.info('Normalising images to mask')

        for i, vol in enumerate(volumes):

            img_a = sitk.GetArrayFromImage(vol)
            mask_a = sitk.GetArrayFromImage(masks[i])

            t = tempfile.TemporaryFile(dir=temp_dir)
            arr_for_mean = np.memmap(t, dtype=img_a.dtype, mode='w+', shape=img_a.shape)

            arr_for_mean[:] = img_a
            arr_for_mean[(mask_a != 1)] = 0
            try:
                # get all values inside mask to calculate mean
                # self.reference_mean = np.mean(img) why is this here anyway
                if fold:
                    # this looks stupid but it stops division by zeroes
                    multi = sitk.MultiplyImageFilter()
                    vol = multi.Execute(vol, self.reference_mean)
                    divis = sitk.DivideImageFilter()
                    volumes[i] = divis.Execute(vol, np.mean(arr_for_mean))
                    #arr = fold_difference * arr  # imagarr = 16bit meandiff = 64bit
                    #tmp = sitk.GetImageFromArray(arr)
                    #tmp.CopyInformation(vol)
                    #volumes[i] = tmp
                else:
                    mean_difference = np.mean(arr_for_mean) - self.reference_mean
                    subtract = sitk.SubtractImageFilter()
                    volumes[i] = subtract.Execute(vol, mean_difference)

            except TypeError:  # Could be caused by imgarr being a short
                # fold difference should not be here
                mean_difference = np.mean(arr_for_mean) - self.reference_mean
                img_a -= int(np.round(mean_difference))
                tmp = sitk.GetImageFromArray(img_a)
                tmp.CopyInformation(vol)
                volumes[i] = tmp


class IntensityHistogramMatch(Normaliser):
    """
    Normalise a set of volumes to the mean of voxel included in a mask.
    In this case each volume needs its mask
    as its not deformed.

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, *kwargs)

    def normalise(self, volumes: List[np.ndarray], ref_vol: np.ndarray):
        """
        Normalises via bin matching to a reference image.
        ThresholdAtMeanIntensityOn() makes

        Parameters
        ----------

        Returns
        -------
        None
            Data is normalised in-place
        """

        logging.info('Using Histogram Matching')
        # logging.info(np.max(ref_vol))

        # ref = sitk.GetImageFromArray(ref_vol)

        # Only need to load the ref volume once

        matcher = sitk.HistogramMatchingImageFilter()
        matcher.SetThresholdAtMeanIntensity(True)
        matcher.SetNumberOfHistogramLevels(256)
        matcher.SetNumberOfMatchPoints(7)
        # matcher.SetReferenceImage(ref_vol)

        for i, img in enumerate(volumes):
            # img = sitk.GetImageFromArray(vol)
            # matcher.SetSourceImage(img)
            volumes[i] = matcher.Execute(img, ref_vol)


class NonRegZNormalise(Normaliser):
    """
    Normalise a set of volumes to the mean of voxel included in a mask.
    In this case each volume needs its mask
    as its not deformed.

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, *kwargs)
        self.reference_mean = None
        self.reference_std = None

    def add_reference(self, ref: np.ndarray, mask: np.ndarray):
        """
        Add the

        Parameters
        ----------

        Returns
        -------
        """
        logging.info('normalising intensity data to mean of the mask')

        # so when we add the reference, we're not storing the image
        # so we can slice it to make computation time quicker
        s = ndimage.find_objects(mask)[0]

        mask = mask[s[0].start:s[0].stop,
               s[1].start:s[1].stop,
               s[2].start:s[2].stop]
        img = ref[s[0].start:s[0].stop,
              s[1].start:s[1].stop,
              s[2].start:s[2].stop]

        # ignore vals outside of mask
        img[(mask != 1)] = 0

        self.reference_mean = np.mean(img).astype(np.uint32)
        self.reference_std = np.std(img).astype(np.uint32)

    def normalise(self, volumes: List[np.ndarray], masks: List[np.ndarray]):
        """
        given paths to registered images, apply linear normalisation so that the mean of the roi across all images are
        the same.

        Create new diretories and place the normalised images in

        Parameters
        ----------

        Returns
        -------
        None
            Data is normalised in-place
        """

        logging.info('Normalising images to mask')

        for i, vol in enumerate(volumes):
            try:
                # get all values inside mask to calculate mean

                img_for_mean = vol

                img_for_mean[(masks[i] != 1)] = 0

                # z-norm from
                # https://www.researchgate.net/publication/229218125_
                # Improving_runoff_prediction_through_the_assimilation_of_the_ASCAT_soil_moisture_product
                vol = (vol - np.mean(img_for_mean).astype(np.uint32)) * (np.std(img_for_mean).astype(np.uint32) /
                                                                         self.reference_std) + self.reference_mean


            except TypeError:  # Could be caused by imgarr being a short
                vol = (vol - np.round(np.mean(img_for_mean))) * (np.round(np.std(img_for_mean))
                                                                 / np.round(self.reference_std)) + np.round(
                    self.reference_mean)


class IntensityMaskNormalise(Normaliser):
    """
    Normalise a set of volumes to the mean of voxe included in a mask.

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, *kwargs)
        self.reference_mean = None

    def add_reference(self, ref: np.ndarray):
        """
        Add the

        Parameters
        ----------

        Returns
        -------
        """
        logging.info('normalising intensity data to mean of the mask')
        self.reference_mean = np.mean(ref)

    def normalise(self, volumes: List[np.ndarray]):
        """
        given paths to registered images, apply linear normalisation so that the mean of the roi across all images are
        the same.

        Create new diretories and place the normalised images in

        Parameters
        ----------

        Returns
        -------
        None
            Data is normalised in-place
        """

        logging.info('Normalising images to mask')

        for vol in volumes:
            try:
                mean_difference = np.mean(vol) - self.reference_mean
                vol -= mean_difference.astype(np.uint16)  # imagarr = 16bit meandiff = 64bit
            except TypeError:  # Could be caused by imgarr being a short
                vol -= int(np.round(mean_difference))


if __name__ == '__main__':
    import argparse

    raise SystemExit('This CLI interafce needs updating')

    parser = argparse.ArgumentParser()
    parser.add_argument('-w', dest='wt', help='wt dir', required=True)
    parser.add_argument('-m', dest='mut', help='mut dir', required=True)
    parser.add_argument('-o', dest='output', help='output dir', required=True)
    parser.add_argument('-s', dest='starts', help='start indices (x, y, z)', required=True, nargs=3, type=int)
    parser.add_argument('-e', dest='ends', help='end indices (x, y, z)', required=True, nargs=3, type=int)

    args = parser.parse_args()
    wt = common.get_file_paths(args.wt)

    mut = common.get_file_paths(args.mut)
    normalise(wt, mut, args.elx_points, args.starts, args.ends)
