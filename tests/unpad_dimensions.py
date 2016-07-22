#!/usr/bin/env python

"""
Jsut a test trying out unpadding images
"""

def unpad(roi_starts, roi_end, config)

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("Get stats on images in a fodler")
    parser.add_argument('-d', '--dir', dest='dir', help='directory with images', required=True)
    parser.add_argument('-o', '--out_file', dest='out_file', help='Optional path to csv file', default=False)
    args = parser.parse_args()
    unpad(args.dir, args.out_file)