#!/usr/bin/env python

import os
import subprocess as sp
import numpy as np
import re
from tempfile import TemporaryFile
import SimpleITK as sitk

minc_dtypes = {'unsigned': {'byte': np.uint8, 'short': np.uint16, 'float': np.float32},
               'signed': {'byte': np.int8, 'short': np.int16, 'float': np.float32}}

def convert(indir, outdir, ext):
    for file_ in os.listdir(indir):
        if not file_.endswith('mnc'):
            continue
        file_ = os.path.abspath(os.path.join(indir, file_))
        info = minc_info(file_)
        mrsg = MincRawSliceGenerator(file_, info)
        vol = mrsg.volume
        base = os.path.splitext(os.path.basename(file_))[0]
        outpath = os.path.join(outdir, base + ext)
        img_out = sitk.GetImageFromArray(vol)
        sitk.WriteImage(img_out, outpath)



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

    def __init__(self, recon, info):
        """The constructor takes a recon path as an argument, and dumps the MINC file to a temporary raw file. The raw file
        is then memory mapped using numpy, from which slices are yielded.

        :param recon: a path to a MINC file.
        """
        super(MincRawSliceGenerator, self).__init__(recon)
        self.ext = 'mnc'

        try:
            self.tmp_file = TemporaryFile(mode='w+b')
            minc_extract = sp.Popen(['mincextract', '-{}'.format(info['dtype']),
                                     '-{}'.format(info['sign']), recon], stdout=self.tmp_file)
            minc_extract.wait()
            self.volume = np.fromfile(self.tmp_file, dtype=info['np_dtype']).reshape(info['shape'])

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

        if 'image:' in line:  # strip non alphanumeric characters

            # Get range
            min_max = re.findall('\d+ to \d+', line)[0]
            info_dict['min'], info_dict['max'] = int(min_max.split()[0]), int(min_max.split()[2])

            regex = re.compile('[^a-zA-Z]')
            info_dict['sign'] = regex.sub('', line.split()[1])
            info_dict['dtype'] = regex.sub('', line.split()[2])

            try:
                info_dict['np_dtype'] = minc_dtypes[info_dict['sign']][info_dict['dtype']]
            except KeyError:
                return None

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

if __name__ == '__main__':
    import sys
    indir = sys.argv[1]
    outdir = sys.argv[2]
    ext = sys.argv[3]
    convert(indir, outdir, ext)