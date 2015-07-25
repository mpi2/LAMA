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


def run(input, output):
    im = sitk.ReadImage(input)
    im_array = sitk.GetArrayFromImage(im)
    shape = im_array.shape
    out_array = np.zeros(shape)

    contrast_weights = _get_contrast_weights([NUMBINS, NUMBINS])

    for z in range(0, shape[0] - (CHUNKSIZE), CHUNKSIZE):
        for y in range(0, shape[1] - (CHUNKSIZE ), CHUNKSIZE):
            for x in range(0, shape[2] - (CHUNKSIZE), CHUNKSIZE):

                chunk = im_array[z: z + CHUNKSIZE, y: y + CHUNKSIZE, x: x + CHUNKSIZE]
                glcm = generate_glcm(chunk)
                #out_array[z: z + CHUNKSIZE, y: y + CHUNKSIZE, x: x + CHUNKSIZE] = _contrast(glcm, contrast_weights)
                out_array[z: z + CHUNKSIZE, y: y + CHUNKSIZE, x: x + CHUNKSIZE] = _stdev(glcm)
    im_out = sitk.GetImageFromArray(out_array)
    sitk.WriteImage(im_out, output)


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
    total = float(np.sum(glcm))
    return glcm/total

def _get_contrast_weights(shape):

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

    """
    return np.std(glcm)



def _contrast(glcm, weights):

    '''Method to create a diagonal matrix'''

    contrast_array = glcm * weights
    contrast = np.sum(contrast_array)
    return contrast


def _digitize(array, bins):
    return bins[np.digitize(array.flatten(), bins) - 1].reshape(array.shape)


if __name__ == '__main__':

    input = '/home/neil/work/glcm_test/20140121RIC8B_15.4_b_wt_rec_scaled_3.1241_pixel_14-1.nrrd'
    output = '/home/neil/work/glcm_test/glcm_var.nrrd'
    array1 = np.random.randint(0, 255, (10, 10, 10))
    array2 = np.random.randint(0, 255, (3, 3))
    run(input, output)




