from lama.utilities import cropper, flipper, lama_convert_16_to_8, lama_pad_volumes
import logging
from pathlib import Path

def run(indir, target_dir):
    # firstly crop volumes
    indir = Path(indir)
    logging.info("Cropping scans")
    cropper.main(indir)

    # cropped files are in a new location
    cropped_dir = indir / "cropped"

    # convert to 8-bit in place (i.e with clobber)
    logging.info("Converting to 8 bit")
    lama_convert_16_to_8.convert_16_bit_to_8bit(cropped_dir,'', clobber=True)

    # flip images - flipper using
    logging.info("flipping")
    flipper.main(cropped_dir)

    # finally pad the volumes with clobber
    logging.info("Padding Volumes")
    lama_pad_volumes.pad_volumes(cropped_dir, '', clobber=True)


def main():
    import argparse

    parser = argparse.ArgumentParser("Arkell Image Processing")
    parser.add_argument('-i', '--input_folder', dest='indirs', help='Raw NRRD directory', required=True,
                        type=str)
    parser.add_argument('-t', '--target_folder', dest='target', help='directory of LAMA target (i.e. pop avg)',
                        required=True,
                        type=str)

    args = parser.parse_args()
    run(args.indirs, args.target)


if __name__ == '__main__':
    main()
