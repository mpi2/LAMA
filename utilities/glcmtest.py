#!/usr/bin/env python

from glcm3d import glcm_run

if __name__ == '__main__':
    import sys
    input = sys.argv[1]
    output = sys.argv[2]
    glcm_run(input, output)