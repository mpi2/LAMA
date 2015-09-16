#!/usr/bin/env python

import sys
import os
import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats





indir = sys.argv[1]
mask = sys.argv[2]


mask_im = sitk.ReadImage(mask)
mask_arr = sitk.GetArrayFromImage(mask_im)



arrays = []
first = True

for impath in os.listdir(indir):

    im = sitk.ReadImage(os.path.join(indir, impath))
    array = sitk.GetArrayFromImage(im)
    #print masked.min, masked.max(0, masked.mean(0))
    print 'hello'
    #array[mask_arr == 0] = np.NAN

    arrays.append(array[mask_arr > 0]._flatten())


stacked = np.vstack(arrays)
stacked = np.sort(stacked, axis=0)


normed = (stacked - stacked.mean(axis=0)).mean(axis=0)


#n = normed.flatten()

plt.plot(normed)
plt.show()

s = stacked.flatten()

nsample = normed.shape[1]
ax4 = plt.subplot(224)
stats.norm.rvs(loc=0, scale=1, size=nsample)
res = stats.probplot(s[~np.isnan(s)], plot=plt)
plt.show()