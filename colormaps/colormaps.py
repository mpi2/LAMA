#!/usr/bin/env python

from __future__ import division
import numpy as np


def generateLUT():
    # Open anatomy file and read organs
    # fAnatomy = open("embryo_anatomy.txt")
    # organs = fAnatomy.readlines()
    # fAnatomy.close()

    # Get last organ number and increment
    #nextNum = int(organs[-2].split(" ")[0]) + 1  # aw yeah one-liner
    nextNum = 0
    # Open output file and write organ data
    fOut = open("tstats_LUT.txt", "w")
    # fOut.writelines(organs[:-1])

    # Interpolate color values for hotred and hotblue
    nextNum = interpolateColors("hotblue.txt", fOut, nextNum, True)
    nextNum = interpolateColors("hotred.txt", fOut, nextNum)

    # Write end of file and close
    # fOut.write(organs[-1])
    fOut.close()


def interpolateColors(colorFile, fOut, nextNum, reverse=0):
    # Get name
    colorName = colorFile.split(".")[0]
    cc = 1;

    # Open color file
    fIn = open(colorFile, "r")
    colors = fIn.readlines()
    fIn.close()

    # Number of colours in each band
    noInterpColors = [25, 25, 25, 25]
    noBaseColors = len(colors) - 1

    # List of colours to write at the end
    IDX = []
    LUT = []

    # Loop through pairs of ARGB colors
    for i in range(noBaseColors):

        argb1 = colors[i].split()
        start = np.array([float(c) for c in argb1])

        argb2 = colors[i + 1].split()
        stop = np.array([float(c) for c in argb2])

        step = (stop - start) / noInterpColors[i]

        for j in range(0, noInterpColors[i]):

            # Compute interpolated colour
            interp = ((start + (j * step)) * 255).astype(int)
            interp = [interp[k] for k in [1, 2, 3, 0]]

            # Convert to string and write to file
            colorString = ' '.join(map(str, interp))
            row = ' {0}_{1} {2}'.format(colorName, cc, colorString)

            # We don't want to include [0 0 0 0] because it's background
            if i + j > 0:
                IDX.append(nextNum)
                LUT.append(row + "\n")
                nextNum += 1
                cc += 1

    # Append the last color because the loop never touches it
    colorString = '255 255 255 255'
    row = ' {0}_{1} {2}'.format(colorName, cc, colorString)
    IDX.append(nextNum)
    LUT.append(row + "\n")
    nextNum += 1

    # Reverse if necessary
    if reverse:
        LUT = reversed(LUT)

    # Write colours
    for idx, color in zip(IDX, LUT):
        fOut.write(str(idx))
        fOut.write(color)


    # Return where we go to
    return nextNum


if __name__ == "__main__":
    generateLUT()