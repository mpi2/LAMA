#!/usr/bin/env python

import sys
import os


dir_ = sys.argv[1]

for subdir, dirs, files in os.walk(dir_):
    for file in files:
        #print os.path.join(subdir, file)
        filepath = subdir + os.sep + file

        if filepath.endswith(".asm"):
            lines = []
            with open(filepath, 'rb') as fh:
                for line in fh:
                    if line.startswith('InitialTransformParametersFileName'):
                        continue
                line.append(line)
            with open(filepath, 'wb') as wh:
                for line in lines:
                    wh.write(line)

