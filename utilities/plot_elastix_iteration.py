#!/usr/bin/env python

import sys

import matplotlib.pyplot as plt

def plot(path):
    data = []

    with open(path) as fh:
        # Get rid of header
        next(fh)
        for line in fh:
            metric = line.split()[1].strip()
            data.append(float(metric))

    plt.plot(data)
    plt.show()


if __name__ == '__main__':
    path = sys.argv[1]
    plot(path)