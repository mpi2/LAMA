#!/usr/bin/env python

import sys

file_ = sys.argv[1]
import matplotlib.pyplot as plt

data = []

with open(file_) as fh:
    for line in fh.readlines():
        data.append(float(line.strip()))

plt.plot(data)
plt.show()


