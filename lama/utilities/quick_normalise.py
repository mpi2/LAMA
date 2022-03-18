

import logging
from pathlib import Path

import numpy as np
from lama import common
import SimpleITK as sitk
import os


def main():
    import argparse
    parser = argparse.ArgumentParser("Run various intensity normalisation methods")
    parser.add_argument('-i', dest='indir', help='directory with vols',
                        required=True)
    parser.add_argument('-l', dest='levels', type=int, help='number of bins within the histogram')
    parser.add_argument('-m', dest='match_num', type=int, help='number of point to match within the histogram')

    parser.add_argument('-o', dest='out_dir', help='output folder for normalised images')



    args = parser.parse_args()

    logging.info("Intensity Normalisation by Histogram bin Matching")
    _dir = Path(args.indir)
    vols = [common.LoadImage(vol) for vol in common.get_file_paths(_dir)]
    names = [os.path.splitext(vol_path.name)[0] for vol_path in common.get_file_paths(_dir)]
    vols.sort()
    names.sort()
    matcher = sitk.HistogramMatchingImageFilter()
    matcher.SetThresholdAtMeanIntensity(True)
    matcher.SetNumberOfHistogramLevels(args.levels)
    matcher.SetNumberOfMatchPoints(args.match_num)

    for i, img in enumerate(vols):
        # img = sitk.GetImageFromArray(vol)
        # matcher.SetSourceImage(img)
        logging.info(f"Normalising {names[i]}")
        vols[i] = matcher.Execute(img, vols[0])
        logging.info(f"Writing Normalised File for {names[i]}")
        file_name = names[i]+".nrrd"
        sitk.WriteImage(vols[i], str(args.outdir / file_name))

if __name__ == '__main__':
    main()