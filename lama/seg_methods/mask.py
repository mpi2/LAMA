#!/usr/bin/env python
import SimpleITK as sitk
from lama.seg_methods.seg import ws_mask, binary_mask, binary_mask_split, auto_thres_mask
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
        parser.add_argument("-external2", dest='external2', action="store_true",
                            help='If external holes are to be filled. New method using 3D whole filling', required=False)
        parser.add_argument("-manual", dest='manual', action="store_true",
                            help='supply manual threshold value', required=False)
        parser.add_argument("-tight", dest='tight', action="store_true",
                            help='If a tight mask is required [Default False]', required=False)
        parser.add_argument("-ndilate", dest='ndilate', action="store_true",
                            help='If no dilation to occur [Default False]', required=False)
        
        # Options for manual binary
        parser.add_argument('-lt', dest='lt', help='[Binary only] level of threshold if man binary used', required=False)
        
        # Options for watershed
        parser.add_argument('-ws_level', dest='ws_level', help='[WS only] watershed level default is 0.9', required=False)
        parser.add_argument('-force', dest='force_background_in', help='[WS only] If the out put for the background '
                                                                       'is not zero', required=False)
        # Options for auto thresh
        parser.add_argument('-bins', dest='bins', help='[auto-thres options only e.g Otsu, Huang] Number of histogram bins'
                                                       ' to use', required=False)
        
        # For option manual binary_split
        parser.add_argument("-split", dest='split', help='If binary_split used. The 3D image is slit in two in its '
                                                         'z-axis by this value', required=False)
        parser.add_argument("-lt1", dest='lt1', help='[Binary_split option only] Lower threshold for first'
                                                     ' half of image', required=False)
        parser.add_argument("-lt2", dest='lt2', help='[Binary_split option only] Lower threshold for second'
                                                     ' half of image', required=False)
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
        
        # This is required if the user wants to perform a binary threshold with a low threshold specified
        if args.lt:
            lt = args.lt
        else:
            lt = False
        
        # This is required if the user wanted to perform a tight mask. (Potentially would get whole volume)
        if args.ws_level:
            ws_level = float(args.ws_level)
        else:
            ws_level = 0.9
        
        # This is for watershed, sometimes need to force the background.
        if not args.force_background_in:
            args.force_background_in = False
        
        # Check if dilation to occur
        if args.ndilate:
            dilate = False
        else:
            dilate = True
        
        #arg.in_file is geetinng a random space from somewhere (which is why I slice it) 

        if args.in_file[-4:] == '.mnc':
            import minc
            sg = minc.MincRawSliceGenerator(args.in_file)
            img = sitk.GetImageFromArray(sg.volume)
        else:
            img = sitk.ReadImage(args.in_file)
        print("I read file")
        
        # All the methods used are in the seg.py module. Split into three types of masking
        
        #====================================================================================================
        # Binary threshold      binary_mask()
        #====================================================================================================
        #   - This includes a binary mask where we choose just a lower threshold value
        #   - Also if no-lower threshold use, we calculate a binary threshold using the maximum inensity and variance of the
        #         the first slice
        
        if option == "binary":
            seg = binary_mask(img,
                              out_dir,
                              lt,
                              tight = args.tight,
                              external = args.external,
                              external2=args.external2,
                              dilate = dilate)
        
        #====================================================================================================
        # Auto threshold        auto_thres_mask()
        #====================================================================================================
        #   - When standard autothres is to be used e.g. otsu, li, huang....
        
        elif option == "otsu" or option == "li" or option == "isodata" or option == "triangle" or option == "huang":
            print("otsu")
            seg = auto_thres_mask(
                option,
                img,
                out_dir,
                bins = bins,
                tight = args.tight,
                external = args.external,
                external2 = args.external2,
                dilate = dilate)

            
        #====================================================================================================
        # Watershed threshold
        #====================================================================================================
        #   - Use the watershed algorithm to mask
        elif option == "ws":
        
            seg = ws_mask(img,
                          out_dir,
                          force_background = args.force_background_in,
                          tight = args.tight,
                          external = args.external,
                          external2=args.external2,
                          ws_level = ws_level,
                          dilate = dilate)
        
        
        #====================================================================================================
        # Binary split threshold
        #====================================================================================================
        #   - Can split the 3D image in the Z-axis by a user defined value
        #   - Can then use a different lower thershold for each section.
        elif option == "binary_split":
        
                seg = binary_mask_split(img,
                                   out_dir,
                                   lt1 = args.lt1,
                                   lt2 = args.lt2,
                                   split = args.split,
                                   tight =args.tight,
                                   external = args.external,
                                   external2=args.external2,
                                   dilate = dilate)
        print("##############################################################################")
        
        
if __name__ == '__main__':
    
    main()
    
    