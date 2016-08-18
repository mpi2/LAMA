#!/usr/bin/env python

import os

folder = '/home/neil/sig/LAMA_results/E18.5/WT_Test_Sets/Test_Set_120216/output/registrations'

find = '(InitialTransformParametersFileName'
replacement = '(InitialTransformParametersFileName "NoInitialTransform")\n'
for dname, dirs, files in os.walk(folder):
    for fname in files:
        if fname == 'TransformParameters.0.txt':
            fpath = os.path.join(dname, fname)
            lines = []
            with open(fpath) as fin:
                for line in fin:
                    if line.startswith(find):
                        line = replacement
                    lines.append(line)

            with open(fpath, "w") as fout:
                for line in lines:
                    fout.write(line)