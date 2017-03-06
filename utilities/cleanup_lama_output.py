#!/usr/bin/env python
"""
Run of a series of lama output directories to cleanup uneeded files
TODO: Everything
"""






if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser("MRC Harwell registration pipeline")
    parser.add_argument('-r', '--root_dir', dest='root_dir', help='Root directory containing lama run', required=True)
    args, _ = parser.parse_known_args()



    args = parser.parse_args()
