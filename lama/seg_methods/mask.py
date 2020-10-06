#!/usr/bin/env python
import SimpleITK as sitk
from lama.seg_methods.seg import auto_thres_mask
import argparse
import textwrap
import os
import sys


#====================================================================================================
# Extracting command line arguments
#====================================================================================================
def main(args):
    # parse arguments using optparse or argparse or what have you
 
        parser = argparse.ArgumentParser(prog='PROG',
                                         formatter_class=argparse.RawDescriptionHelpFormatter,
                                         description='''Masking Script ''',
                                         epilog=textwrap.dedent('''\
                                         -------------------------------------------------------------------------
        
                                         Example Usage:
        
                                         -Default mask run.
                                         python mask.py -i [infile]
        
                                         -Example mask with otsu, external hole filling and a tight mask
                                         python mask.py -i [infile] -o [out folder] -option otsu -external -tight
        
                                         -Example mask with huang, external hole filling and rough mask (default)
                                         python mask.py -i [infile] -o [out folder] -option huang -external
        
                                         -Example mask with ws, external hole filling and rough mask (default)
                                         python mask.py -i [infile] -o [out folder] -option ws
                                         '''))
        
        # Global options
        parser.add_argument('-i', dest='in_file', help='uCT recon image', required=True)
        parser.add_argument('-o', dest='out_dir', help='destination for output', required=False)
        parser.add_argument('-option', dest='option', help='Otsu, huang, li, ws, binary, binary_split', required=False)
        parser.add_argument("-external", dest='external', action="store_true",
                            help='If external holes are to be filled.', required=False)


        parser.add_argument("-tight", dest='tight', action="store_true",
                            help='If a tight mask is required [Default False]', required=False)
        parser.add_argument("-ndilate", dest='ndilate', action="store_true",
                            help='If no dilation to occur [Default False]', required=False)
        

        # Options for auto thresh
        parser.add_argument('-bins', dest='bins', help='[auto-thres options only e.g Otsu, Huang] Number of histogram bins'
                                                       ' to use', required=False)
        args = parser.parse_args(args)
        #====================================================================================================
        # Setting up folders and parameters
        #====================================================================================================
        print("################################ Masking Script ################################")
        
        print("creating folders")
        # Create a folder if one is not already been given by the user
        if args.out_dir:
            out_dir = args.out_dir
        else:
            basename = os.path.basename(args.in_file)
            path, file = os.path.split(args.in_file)
            out_dir = os.path.join(path,os.path.splitext(basename)[0])
        
        if not os.path.exists(out_dir):
                os.makedirs(out_dir)
        
        # Get the option of mask from the user. Defaults to binary
        
        if args.option:
            option = args.option
        else:
            option = "binary"
        
        # Get the option of mask from the user. Defaults to binary
        if args.bins:
            bins = args.bins
        else:
            bins = False

        
        # Check if dilation to occur
        if args.ndilate:
            dilate = False
        else:
            dilate = True
        
        #arg.in_file is geetinng a random space from somewhere (which is why I slice it) 
        img = sitk.ReadImage(args.in_file)

        #====================================================================================================
        # Auto threshold        auto_thres_mask()
        #====================================================================================================
        #   - When standard autothres is to be used e.g. otsu, li, huang....
        
        if option == "otsu" or option == "li" or option == "isodata" or option == "triangle" or option == "huang":
            seg = auto_thres_mask(
                option,
                img,
                out_dir,
                bins = bins,
                tight = args.tight,
                external = args.external,
                external2 = args.external2,
                dilate = dilate)

if __name__ == '__main__':
    main()
