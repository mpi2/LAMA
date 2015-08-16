#!/usr/bin/python

"""
Implement a 3D GLCM for use in the registration  pipeline

"""


import numpy as np
import SimpleITK as sitk
from multiprocessing import Process, Queue
import tempfile
MAXINTENSITY = 255
import copy
from os.path import join
import uuid



class GlcmGenerator(Process):
    """
    Currently only works on 8bit images
    """
    def __init__(self, task_queue, result_queue, chunksize, numbins, mask=None):
        super(GlcmGenerator, self).__init__()
        self.result_queue = result_queue

        self.task_queue = task_queue

        #img_shape = sitk.GetArrayFromImage(sitk.ReadImage(img_path)).shape

        if mask:
            self.mask = sitk.GetArrayFromImage(sitk.ReadImage(mask))
        #     if self.mask.shape != img_shape:
        #         raise ValueError("Mask and input image shape need to be the same")
        # else:
            self.mask = None

        self.chunksize = chunksize
        self.numbins = numbins
        self.bins = np.array(range(0, 256, 256/numbins))

    def run(self):
        """
        Multithread this stage as it's the most time consuming
        :return:
        """
        while True:
            img_path = self.task_queue.get()
            if img_path is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break
            self.task_queue.task_done()
            tmp_dir = tempfile.gettempdir()
            prefix = str(uuid.uuid4())
            tmp_file = join(tmp_dir, prefix)
            img_glcms = self._generate_glcms(img_path)
            np.save(tmp_file, img_glcms)
            tempfile_ext = join(tmp_file + '.npy')
            print 'gen_glcm', prefix
            self.result_queue.put(tempfile_ext)
        return

    def _generate_glcms(self, im_path):
        """
        Generate 3D chunks from an image array and generate a glcm for each chunk. If a mask is avaiable set the
        glcm at fully masked regions  to 'None'
        """
        im = sitk.ReadImage(im_path)
        im_array = sitk.GetArrayFromImage(im)
        shape = im_array.shape
        img_glcms = []

        for z in range(0, shape[0] - self.chunksize, self.chunksize):
            for y in range(0, shape[1] - self.chunksize, self.chunksize):
                for x in range(0, shape[2] - self.chunksize, self.chunksize):

                    chunk = im_array[z: z + self.chunksize, y: y + self.chunksize, x: x + self.chunksize]
                    if self._outside_of_mask(z, y, x):
                        glcm = _generate_glcm(chunk, self.numbins, self.bins)
                    else:
                        glcm = None
                    img_glcms.append(glcm)
        return img_glcms

    def _outside_of_mask(self, z, y, x):
        """
        Check whether an roi is part of the mask. If masdk not set, return True

        Parameters
        ----------
        z, y, x: int
            the start indices of the roi
        """
        if self.mask != None:
            roi = self.mask[z: z + self.chunksize, y: y + self.chunksize, x: x + self.chunksize]
            if np.any(roi):  # Any values other than zero
                return True
            else:
                return False
        else:
            return True


def _generate_glcm(array, num_bins, bins):
    """
    Currently just using a pixel one away x and y. Try the invariant direction. ie the 6 neighbouring pixels
    It will be slow, so probably do this as a c extension
    :param array:
    :return:
    """
    glcm = np.zeros((num_bins, num_bins))

    it = np.nditer(array, flags=['multi_index'])
    while not it.finished:
        reference = array[it.multi_index]
        try:
            neighbour1 = array[it.multi_index[0], it.multi_index[1] + 2, it.multi_index[2] + 2]
        except IndexError:
            pass  # At the edge of the array
        try:
            neighbour2 = array[it.multi_index[0], it.multi_index[1] - 2, it.multi_index[2] - 2]
        except IndexError:
            pass # At the edge of the array
        else:
            co_pair1 = np.digitize([reference, neighbour1], bins) - 1
            co_pair2 = np.digitize([reference, neighbour2], bins) - 1
            try:
                glcm[co_pair1[0], co_pair1[1]] += 1
                glcm[co_pair1[1], co_pair1[0]] += 1  # We need a symetrical array
                glcm[co_pair2[0], co_pair2[1]] += 1
                glcm[co_pair2[1], co_pair2[0]] += 1  # We need a symetrical array
            except IndexError:
                pass
        it.iternext()

    norm_glcm = _normalize(glcm)
    return norm_glcm


def _normalize(glcm):
    """
    Replace glcm counts with probabilities
    """
    total = float(np.sum(glcm))
    return glcm/total


def _stdev(glcm):
    """
    Gets standard deviation metric of a glcm
    """
    return np.std(glcm)


def _discretize(array, bins):
    """
    bin array intensity values
    """
    return bins[np.digitize(array, bins) - 1]


class TextureCalculations(object):
    def __init__(self, glcms, glcm_info):
        """
        Parameters
        ----------
        glcms: ndarray
            ijk array. i: number of specimens, j: number of glcms, k: shape of glcm
            k can also be 'None' if for example it's a msked region
        glcm_info: dict
            info about the glcm. chunksize etc
        """
        self.glcms = glcms
        self.glcm_info = glcm_info

    def _get_texture(self):
        raise NotImplementedError

    def get_results(self):
        """
        Return the texture results.
        """
        results = []
        for glcm in self.glcms:
            # If the region was outside of the mask, the GLCM for this region will be 'None'
            if glcm != None:
                results.append(self._get_texture(glcm))
            else:
                results.append(0)
        return results


class ContrastTexture(TextureCalculations):
    def __init__(self, *args):
        super(ContrastTexture, self).__init__(*args)
        numbins = self.glcm_info['num_bins']
        self.contrast_weights = self._get_contrast_weights((numbins, numbins))

    def _get_texture(self, glcm):
        """
        Get contrast measure of a glcm
        """
        contrast_array = glcm * self.contrast_weights
        contrast = np.sum(contrast_array)
        return contrast

    @staticmethod
    def _get_contrast_weights(shape):
        """
        Get the contrast weight metrix. Increased weight with increased difference between x and y
        """

        size = shape[0] * shape[1]
        weights = np.zeros(size).reshape(shape)
        it = np.nditer(weights, flags=['multi_index'])
        while not it.finished:
            index = it.multi_index
            diff = abs(index[0] - index[1])
            w = diff * diff
            weights[index[0], index[1]] = w
            it.iternext()
        return weights


class ASMTexture(TextureCalculations):
    """
    http://www.fp.ucalgary.ca/mhallbey/asm.htm
    """
    def __init__(self, *args):
        super(ASMTexture, self).__init__(*args)

    def _get_texture(self, glcm):
        if glcm != None:
            asm = np.square(glcm).sum()
        else:
            asm = None
        return asm

class EntropyTexture(TextureCalculations):
    """
    """
    def __init__(self, *args):
        super(EntropyTexture, self).__init__(*args)

    def _get_texture(self, glcm):
        """
        Not sure if the entropy calculation is correct.
        For arrays len(100) with same value, I get an entropy of ~1
        For same array with random values I get entropy of >2
        """
        e = glcm.flatten()
        entropy = np.sum([p*np.log2(1.0/p) for p in e])
        return entropy




if __name__ == '__main__':
    import sys
    input_ = sys.argv[1]
    output = sys.argv[2]

    g = GlcmGenerator(input_)
    g.write_contrast_glcm_image(output)




