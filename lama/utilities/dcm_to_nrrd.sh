#!/bin/bash

#DICOM to NRRD  file  bash script

#TODO: add genotype detection by identifying 'wt' in the parent directory 
#(uCT scans are performed before genotyping to reduce bias and therefore 
#a method of simplified labelling post scanning will be required)
#May be able to intergrate code from tiff_to_minc.sh

#loops through folders
mkdir nrrd_out

for directory in */;
do 
    #Go within the specific DICOM directory:
    dir_name=${directory%/*// /_}
    cd ${dir_name}

    #Error trap for spaces in dicom filenames from previous PhD student
    for f in *\ *; do mv "$f" "${f// /_}"; done

    #Make directory and perform the conversion
    cd ../
    /home/minc/.config/NA-MIC/Extensions-29402/SlicerDcm2nii/lib/Slicer-4.11/qt-scripted-modules/Resources/bin/dcm2niix -1 -d 0 -l y -f ${dir_name%/} -o nrrd_out -e y -z n ${dir_name}

    #Do some basic clean-up - we don't care about .json files
    rm nrrd_out/*.json
done




















