#!/bin/bash

#Bash script to convert DICOM files to NRRD format in batches. 
#Currenty converts DICOMs to 16-bit utype (and as such may need modifiying).  

#Dependencies: 
#     Slicer
#     dcm2niix module for Slicer (path to executable module depends on the host) 


#TODO: add genotype detection by identifying 'wt' in the parent directory 
#(uCT scans are performed before genotyping to reduce bias and therefore 
#a method of simplified labelling post scanning will be required)

#TODO: double check headers are compatible with LAMA. 

#Make directory
mkdir nrrd_out

#loops through folders
for directory in */;
do 
    #Go within the specific DICOM directory:
    dir_name=${directory%/*// /_}
    cd ${dir_name}

    #Error trap for spaces in dicom filenames from previous PhD student
    for f in *\ *; do mv "$f" "${f// /_}"; done

    #Perform the conversion
    #TODO: find code to identify where dcm2niix is located. 
    cd ../
    /home/minc/.config/NA-MIC/Extensions-28257/SlicerDcm2nii/lib/Slicer-4.10/qt-scripted-modules/Resources/bin/dcm2niix -1 -d 0 -f ${dir_name%/} -o nrrd_out -e y -z n ${dir_name} 
done




















