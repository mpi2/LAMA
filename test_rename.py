#!/usr/bin/python3

import os
import shutil

def get_file_paths(folder, extension_tuple=('.nrrd', '.tiff', '.tif', '.nii', '.bmp', 'jpg', 'mnc', 'vtk', 'bin'),
                   pattern=None, ignore_folder=""):
    """
    Test whether input is a folder or a file. If a file or list, return it.
    If a dir, return all images within that directory.
    Optionally test for a pattern to sarch for in the filenames
    """

    if not os.path.isdir(folder):
        return False
    else:
        paths = []
        for root, subfolders, files in os.walk(folder):
            if ignore_folder in subfolders:
                subfolders.remove(ignore_folder)
            for filename in files:
                if filename.lower().endswith(extension_tuple):
                    if pattern:
                        if pattern and pattern not in filename:
                            continue
                    #paths.append(os.path.abspath(os.path.join(root, filename))) #Broken on shared drive
                    paths.append(os.path.join(root, filename))
        return paths
indir = '/home/neil/bit/LAMA_results/E14.5/paper_runs/output/inverted_labels'

file_paths = get_file_paths(indir)
for p in file_paths:
	new = p.replace('seg_', '')
	os.rename(p, new)