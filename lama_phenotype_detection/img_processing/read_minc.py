# Copyright 2016 Medical Research Council Harwell.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#
# @author James Brown <james.brown@har.mrc.ac.uk>

import subprocess as sp
import numpy as np
import re
from tempfile import NamedTemporaryFile

minc_dtypes = {'unsigned': {'byte': np.uint8, 'short': np.uint16, 'float': np.float32},
               'signed': {'byte': np.int8, 'short': np.int16, 'float': np.float32}}

def minc_to_numpy(minc_file):

    info = minc_info(minc_file)
    if not info:
        return False
    mrsg = MincRawSliceGenerator(minc_file, info)
    vol = mrsg.volume
    return vol


def mincstats_to_numpy(minc_file):

    info = minc_info(minc_file)
    if not info:
        return False
    info['dtype'] = 'float'  # have to force it, because mincinfo is hopeless
    info['np_dtype'] = np.float32
    mrsg = MincRawSliceGenerator(minc_file, info)
    vol = mrsg.volume
    return vol


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

        tmp_file = NamedTemporaryFile() # TemporaryFile() seems not to work with Python3.4

        sp.call(['mincextract', '-{}'.format(info['dtype']),
                             '-{}'.format(info['sign']), recon], stdout=tmp_file)

        self.volume = np.fromfile(tmp_file.name, dtype=info['np_dtype']).reshape(info['shape'])

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

    try:
        info = sp.check_output(['mincinfo', recon], universal_newlines=True)
    except OSError as e:
        raise OSError("Minc tools not installed\n{}".format(e))
    #info = str(info)

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
