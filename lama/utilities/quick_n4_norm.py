

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
    #parser.add_argument('-l', dest='levels', type=int, help='number of bins within the histogram')
    #parser.add_argument('-m', dest='match_num', type=int, help='number of point to match within the histogram')

    parser.add_argument('-o', dest='out_dir', help='output folder for normalised images')

    #parser.add_argument('-r', dest='ref_vol', help='Path of Reference Volume To Use')



    args = parser.parse_args()

    logging.info("Intensity Normalisation by N4")

    print(Path(args.indir))

    _files = common.get_file_paths(Path(args.indir))

    _files.sort(key=lambda x: os.path.basename(x))

    # should be labelling properly now
    vols = [common.LoadImage(_path).img for _path in _files]

    names = [os.path.splitext(vol_path.name)[0] for vol_path in _files]


    Otsu = sitk.OtsuThresholdImageFilter()
    dilate = sitk.BinaryDilateImageFilter()
    downsampler = sitk.ShrinkImageFilter()
    N4 = sitk.N4BiasFieldCorrectionImageFilter()

    #if args.ref_vol:
    #    ref_vol = common.LoadImage(Path(args.ref_vol)).img
    #else:
    #    ref_vol = vols[0]

    for i, vol in enumerate(vols):
        print(i)
        logging.info("Getting mask using via the otsu algorithm")
        inv_mask = Otsu.Execute(vol)
        o_mask = sitk.InvertIntensity(inv_mask, 1)

        o_mask = sitk.ConnectedComponent(o_mask != o_mask[0, 0, 0])

        # sitk.WriteImage(seg, os.path.join(output, name + "_all_connected.nrrd"))
        o_mask = sitk.RelabelComponent(o_mask)
        o_mask = o_mask == 1
        # sitk.WriteImage(seg, os.path.join(output, name + "_largest_connected.nrrd"))

        # lets see if dilate with a tight kernal fixes getting stupid dots everywhere.

        dilate.SetKernelRadius([1, 1, 1])
        dilate.SetKernelType(sitk.sitkBall)
        o_mask = dilate.Execute(o_mask)
        o_mask.CopyInformation(vol)

        logging.info('Using N4 bias correction')

        # downsample images

        down_sampled_img = downsampler.Execute(vol)

        down_sampled_mask = downsampler.Execute(o_mask)

        N4_vol = N4.Execute(down_sampled_img, down_sampled_mask)
        log_bias_field = N4.GetLogBiasFieldAsImage(vol)
        vols[i] = vol / sitk.Exp(log_bias_field)
        logging.info(f"Writing Normalised File for {names[i]}")

        file_name = names[i]+".nrrd"
        sitk.WriteImage(vols[i], str(Path(args.out_dir)/ file_name))

if __name__ == '__main__':
    main()