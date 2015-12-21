#! /usr/bin/env python

"""

"""

import SimpleITK as sitk
from SimpleITK import ReadImage, WriteImage, GetArrayFromImage, GetImageFromArray
import numpy as np
from os.path import join
import sys
import matplotlib.pyplot as plt
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common


def run(lama_data_path, cutoff):
    # im = ReadImage(lama_data_path)
    # arr = GetArrayFromImage(im)

    # Try out gettiong features along z axis
    data = np.load(lama_data_path)
    qvals = data['qvals'][0]
    tvals = data['tvals'][0]

    #filter the tstats
    tvals[qvals >= cutoff] = 0

    neg_features = np.copy(tvals)
    neg_features[neg_features > 0] = 0
    pos_features = np.copy(tvals)
    pos_features[pos_features < 0] = 0



    plt.plot(list(np.mean(neg_features, axis=(1,2))))
    plt.title("Negative values")
    plt.show()

    plt.plot(list(np.mean(pos_features, axis=(1,2))))
    plt.title("Positive values")
    plt.show()




if __name__ == '__main__':
    lama_data_path = sys.argv[1]
    cutoff = 0.05
    run(lama_data_path, cutoff)
