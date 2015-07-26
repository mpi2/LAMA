#!/usr/bin/python



"""
Implement a 3D GLCM for use in the registration  pipeline




"""


import numpy as np
import SimpleITK as sitk

NUMBINS = 16
BINS = np.array(range(0, 256, NUMBINS))
MAXINTENSITY = 255
CHUNKSIZE = 10

class Glcm(object):
    """
    Currently only works on 8bit images
    """
    def __init__(self, img_in, numbins=16, threads=4):
        self.img_in = img_in
        im = sitk.ReadImage(img_in)
        self.im_array = sitk.GetArrayFromImage(im)
        self.numbins = numbins
        self.threads = threads
        self.img_glcms = self.generate_glcms(self.im_array)

    def contrast_glcm(self, img_out):
        print 'getting contrast'
        contrast_weights = _get_contrast_weights((self.numbins, self.numbins))
        contrasts = []
        for glcm in self.img_glcms:
            contrasts.append(_contrast(glcm, contrast_weights))
        self.write_results(contrasts, img_out)

    def write_results(self, result, path):
        shape = self.im_array.shape
        out_array = np.zeros(shape)

        i = 0
        
        for z in range(0, shape[0] - CHUNKSIZE, CHUNKSIZE):
            print 'w', z
            for y in range(0, shape[1] - CHUNKSIZE, CHUNKSIZE):
                for x in range(0, shape[2] - CHUNKSIZE, CHUNKSIZE):
                    out_array[z: z + CHUNKSIZE, y: y + CHUNKSIZE, x: x + CHUNKSIZE] = result[i]
                    i += 1
        out_img = sitk.GetImageFromArray(out_array)
        sitk.WriteImage(out_img, path)

    def generate_glcms(self, im_array):

        shape = im_array.shape
        img_glcms = []

        #contrast_weights = _get_contrast_weights([NUMBINS, NUMBINS])

        for z in range(0, shape[0] - (CHUNKSIZE), CHUNKSIZE):
            print z
            for y in range(0, shape[1] - (CHUNKSIZE ), CHUNKSIZE):
                for x in range(0, shape[2] - (CHUNKSIZE), CHUNKSIZE):

                    chunk = im_array[z: z + CHUNKSIZE, y: y + CHUNKSIZE, x: x + CHUNKSIZE]
                    glcm = generate_glcm(chunk)
                    #out_array[z: z + CHUNKSIZE, y: y + CHUNKSIZE, x: x + CHUNKSIZE] = _contrast(glcm, contrast_weights)
                    #img_glcms[z: z + CHUNKSIZE, y: y + CHUNKSIZE, x: x + CHUNKSIZE] = _stdev(glcm)
                    img_glcms.append(glcm)
        return img_glcms
        #im_out = sitk.GetImageFromArray(img_glcms)
        #sitk.WriteImage(im_out, output_img)


def generate_glcm(array):
    """
    Currently just using a pixel one away x and y. Try the invariant direction. ie the 6 neighbouring pixels
    It will be slow, so probably do this as a c extension
    :param array:
    :return:
    """
    glcm = np.zeros((NUMBINS, NUMBINS))
    # binned = _digitize(array, BINS)

    it = np.nditer(array, flags=['multi_index'])
    while not it.finished:
        try:
            reference = array[it.multi_index]
            neighbour = array[it.multi_index[0], it.multi_index[1] + 1, it.multi_index[2] + 1]
        except IndexError:
            pass # At the edge of the array
        else:
            co_pair = np.digitize([reference, neighbour], BINS) - 1
            glcm[co_pair[0], co_pair[1]] += 1
            glcm[co_pair[1], co_pair[0]] += 1 # We need a symetrical array
        it.iternext()

    norm_glcm = _normalize(glcm)
    return norm_glcm


def _normalize(glcm):
    """
    Replace glcm counts with probabilities
    """
    total = float(np.sum(glcm))
    return glcm/total

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

def _stdev(glcm):
    """
    Gets standard deviation metric of a glcm
    """
    return np.std(glcm)



def _contrast(glcm, weights):
    """
    Get contrast measure of a glcm
    """

    contrast_array = glcm * weights
    contrast = np.sum(contrast_array)
    return contrast


def _digitize(array, bins):
    """
    bin array intensity values
    """
    return bins[np.digitize(array.flatten(), bins) - 1].reshape(array.shape)


if __name__ == '__main__':
    import sys
    input_ = sys.argv[1]
    output = sys.argv[2]

    g = Glcm(input_)
    g.contrast_glcm(output)




