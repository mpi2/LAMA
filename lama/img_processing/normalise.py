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
            elif type_ == 'N4biascorrection':
                return IntensityN4Normalise()
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

    def get_all_wt_vols_and_masks(self, _dir):
        baseline_dir = Path(_dir).parent / "baseline"
        vol_paths = [path for path in common.get_file_paths(baseline_dir) if "rigid" in str(path)]
        mask_paths = [path for path in common.get_file_paths(baseline_dir) if "inverted_stats_mask" in str(path)]

        vol_paths.sort(key=lambda x: os.path.basename(x))
        mask_paths.sort(key=lambda x: os.path.basename(x))

        vols = [common.LoadImage(_path) for _path in vol_paths]
        masks = [common.LoadImage(_path) for _path in mask_paths]

        return vols, masks



    def gen_otsu_masks(self, volumes: List[np.ndarray], file_names: List[Path]=None):
        '''
        Creates an otsu for each scan
        Parameters
        ----------
        volumes - list of volumes

        Returns
        -------

        '''
        logging.info("Creating_otsu_masks")
        o_masks = [None]* len(volumes)
        print(len(volumes))
        if ~isinstance(volumes, list):
            # stops code from breaking in radiomics runner
            volumes = [volumes]
        for i, vol in enumerate(volumes):
            print(vol)
            Otsu = sitk.OtsuThresholdImageFilter()

            inv_mask = Otsu.Execute(vol)
            o_mask = sitk.InvertIntensity(inv_mask, 1)

            o_mask = sitk.ConnectedComponent(o_mask != o_mask[0, 0, 0])

            # sitk.WriteImage(seg, os.path.join(output, name + "_all_connected.nrrd"))
            o_mask = sitk.RelabelComponent(o_mask)
            o_mask = o_mask == 1
            # sitk.WriteImage(seg, os.path.join(output, name + "_largest_connected.nrrd"))

            # lets see if dilate with a tight kernal fixes getting stupid dots everywhere.
            dilate = sitk.BinaryDilateImageFilter()
            dilate.SetKernelRadius([1, 1, 1])
            dilate.SetKernelType(sitk.sitkBall)
            o_masks[i] = dilate.Execute(o_mask)
            o_masks[i].CopyInformation(vol)

            #o_dir = file_names[0].parent.parent / "otsu_thresholds"
            #os.makedirs(o_dir, exist_ok=True)
            #sitk.WriteImage(o_mask, str(Path(o_dir) / os.path.basename(file_names[i])))
        return o_masks

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
        means = []
        for i, vol in enumerate(ref):
            img = sitk.GetArrayFromImage(ref[i].img)
            mask = sitk.GetArrayFromImage(ref_mask[i].img)

            s = ndimage.find_objects(mask)[0]

            mask = mask[s[0].start:s[0].stop,
                   s[1].start:s[1].stop,
                   s[2].start:s[2].stop]
            img = img[s[0].start:s[0].stop,
                  s[1].start:s[1].stop,
                  s[2].start:s[2].stop]

            # test if this improves speed

            # ignore vals outside of mask
            img = img[mask == 1]

            means.append(np.mean(img))


        self.reference_mean = np.mean(means)

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
            img_a = sitk.GetArrayFromImage(vol.img)
            mask_a = sitk.GetArrayFromImage(masks[i].img)
            t = tempfile.TemporaryFile(dir=temp_dir)
            img_a = img_a[mask_a == 1]
            arr_for_mean = np.memmap(t, dtype=img_a.dtype, mode='w+', shape=img_a.shape)
            arr_for_mean[:] = img_a

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
                    img = vol.img
                    volumes[i] = subtract.Execute(img, float(mean_difference))

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

        #try:
        #    ref_vol_path = Path(config.config_dir / config['reference_vol'])
        #    self.ref_vol = common.LoadImage(ref_vol_path)
        #except KeyError:
        #    self.ref_vol = None

    def normalise(self, volumes: List[np.ndarray], ref_vol: np.ndarray = None):
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

        # Get the Population average as the ref vol if not provided.
        #ref_vol = self.ref_vol if self.ref_vol else ref_vol

        # Only need to load the ref volume once
        matcher = sitk.HistogramMatchingImageFilter()
        matcher.SetThresholdAtMeanIntensity(True)

        for i, img in enumerate(volumes):
            try:
                volumes[i] = matcher.Execute(img, ref_vol)
            except RuntimeError: # needs casting
                #img = sitk.Cast(img, sitk.sitkFloat32)
                #ref_vol = sitk.Cast(ref_vol, sitk.sitkFloat32)
                volumes[i] = matcher.Execute(img, ref_vol)

class IntensityN4Normalise(Normaliser):
    """
    Use N4 normalisation to normalise images - needs to have a mask specified or else all values > 1
    are masked.

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, *kwargs)

        #try:
        #    ref_vol_path = Path(config.config_dir / config['reference_vol'])
        #    self.ref_vol = common.LoadImage(ref_vol_path)
        #except KeyError:
        #    self.ref_vol = None
    def gen_otsu_masks(self, vol: List[np.ndarray], file_names: List[Path]=None):
        '''
        Creates an otsu for each scan
        Parameters
        ----------
        volumes - list of volumes

        Returns
        -------

        '''
        logging.info("Creating_otsu_masks")

        Otsu = sitk.OtsuThresholdImageFilter()

        inv_mask = Otsu.Execute(vol)
        o_mask = sitk.InvertIntensity(inv_mask, 1)

        o_mask = sitk.ConnectedComponent(o_mask != o_mask[0, 0, 0])

        # sitk.WriteImage(seg, os.path.join(output, name + "_all_connected.nrrd"))
        o_mask = sitk.RelabelComponent(o_mask)
        o_mask = o_mask == 1
        # sitk.WriteImage(seg, os.path.join(output, name + "_largest_connected.nrrd"))

        # lets see if dilate with a tight kernal fixes getting stupid dots everywhere.
        dilate = sitk.BinaryDilateImageFilter()
        dilate.SetKernelRadius([1, 1, 1])
        dilate.SetKernelType(sitk.sitkBall)
        o_mask = dilate.Execute(o_mask)
        o_mask.CopyInformation(vol)

            #o_dir = file_names[0].parent.parent / "otsu_thresholds"
            #os.makedirs(o_dir, exist_ok=True)
            #sitk.WriteImage(o_mask, str(Path(o_dir) / os.path.basename(file_names[i])))
        return o_mask

    def normalise(self, img: List[np.ndarray], mask=List[np.ndarray]):
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

        logging.info('Using N4 bias correction')

        # downsample images
        downsampler = sitk.ShrinkImageFilter()

        down_sampled_img = downsampler.Execute(img)

        down_sampled_mask = downsampler.Execute(mask)

        N4 = sitk.N4BiasFieldCorrectionImageFilter()

        N4_vol = N4.Execute(down_sampled_img, down_sampled_mask)

        log_bias_field = N4.GetLogBiasFieldAsImage(img)
        img = img / sitk.Exp(log_bias_field)
        sitk.WriteImage(img, "E:/220607_two_way/radiomics_output/test_b4.nrrd")



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
