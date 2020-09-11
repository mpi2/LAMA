__author__ = 'james'

import numpy as np
import subprocess as sp
from tempfile import TemporaryFile
import re

class SliceGenerator(object):

    def __init__(self, recon):
        self.recon = recon
        self.slice_index = 0

    def slices(self):
        """The slices method should yield xy image slices from a memory mapped numpy array."""
        raise NotImplementedError("Ths method needs overriding")

    def dtype(self):
        """The dtype method should return the datatype of the memory mapped numpy array"""
        raise NotImplementedError("Ths method needs overriding")

    def shape(self):
        """The shape method should return the shape of the memory mapped numpy array in x, y, z order."""
        raise NotImplementedError("Ths method needs overriding")

class MincRawSliceGenerator(SliceGenerator):
    """The MincRawSliceGenerator class extends SliceGenerator, yielding slices from a single MINC (Medical Image NetCDF)
    file, having been dumped to a temporary raw file via minctoraw. mincinfo is used to determine the file type/dimensions
    """

    def __init__(self, recon):
        """The constructor takes a recon path as an argument, and dumps the MINC file to a temporary raw file. The raw file
        is then memory mapped using numpy, from which slices are yielded.

        :param recon: a path to a MINC file.
        """
        super(MincRawSliceGenerator, self).__init__(recon)
        self.ext = 'mnc'

        try:
            self.tmp_file = TemporaryFile(mode='w+b')
            info = minc_info(recon)
            minc_extract = sp.Popen(['mincextract', '-{}'.format(info['dtype']),
                                     '-{}'.format(info['sign']), recon], stdout=self.tmp_file)
            minc_extract.wait()
            self.volume = np.memmap(self.tmp_file, dtype=np.uint16, shape=info['shape'], mode='r')

        except Exception as e:
            print e
            raise ReconFormatError("Error opening file as MINC")

    def slices(self, start=0):
        """Slices are yielded one slice at a time from the memory mapped numpy array
        """
        try:
            for i in range(self.volume.shape[0]):
                yield self.volume[i, :, :]
        except Exception:
            raise CorruptReconError("Error yielding slices from MINC file")

    def dtype(self):
        """Overrides the superclass to return the data type of the MINC file i.e. 8 bit/16 bit.
        """
        return self.volume.dtype

    def shape(self):
        """Overrides the superclass to return the shape of the MINC file.
        """
        return self.volume.shape[::-1]


def minc_info(recon):

    info = sp.check_output(['mincinfo', recon])

    info_dict = {}
    dims = []
    for line in info.splitlines():

        if 'image:' in line:
            info_dict['sign'] = re.sub('[^0-9a-zA-Z]+', '', line.split()[1])
            info_dict['dtype'] = re.sub('[^0-9a-zA-Z]+', '', line.split()[2])
        elif 'dimensions' not in line and any(space in line for space in ['xspace', 'yspace', 'zspace']):
            spacing = line.split()
            dims.append(int(spacing[1]))
            info_dict['voxel_size'] = float(spacing[2]) * 1000  # in microns

    info_dict['shape'] = tuple(dims)  # zspace, yspace, xspace
    return info_dict

class ReconFormatError(Exception):

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


class CorruptReconError(Exception):

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)
