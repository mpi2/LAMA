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

from lama.paths import specimen_iterator
from lama import common
from lama.registration_pipeline.validate_config import LamaConfig

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
            if isinstance(type_ ,(list, )):  # Not working at the moment
                if len(type_) != 3:
                    return None
                return None # RoiNormalise

            elif type_ == 'mask':
                return IntensityMaskNormalise()

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
            config_file = common.getfile_endswith('.toml') # Get the Lama config from the specimen directory
            config = LamaConfig(config_file )
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
