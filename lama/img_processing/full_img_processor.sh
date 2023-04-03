#!/bin/bash

#DICOM to NRRD  file  bash script

#TODO: add genotype detection by identifying 'wt' in the parent directory 
#(uCT scans are performed before genotyping to reduce bias and therefore 
#a method of simplified labelling post scanning will be required)
#May be able to intergrate code from tiff_to_minc.sh

# get the latest dcm to niix and unzip - ToDO: decide whether to just download the zip or proper install... 

curl -fLO https://github.com/rordenlab/dcm2niix/releases/latest/download/dcm2niix_lnx.zip
unzip dcm2niix_lnx.zip

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

    #TODO: check if -l o/n is better!!!!
    ./dcm2niix -1 -d 0 -f ${dir_name%/} -o nrrd_out -e y -z n ${dir_name}

    #Do some basic clean-up - we don't care about .json files
    rm nrrd_out/*.json
done


# img processor.sh
mv cropper.py nrrd_out/


cd nrrd_out

mkdir cropped masked

python3 cropper.py

# 16-bit to 8-bit conversion

mkdir converted

lama_convert_16_to_8 -i cropped -o converted

# flipping scans to matc the population average

mv ../flipper.py converted/

cd converted

python3 flipper.py

# padding
cd ../

mkdir needs_padding
mkdir padded
cp -r ../target needs_padding/

mv converted needs_padding/

lama_pad_volumes -i needs_padding -o padded

mv padded ../


cd padded

# make folders
for cond in baseline, mutants, treatment, mut_treatment;
do
	mkdir ${cond}
	mkdir ${cond}/inputs
	mkdir ${cond}/inputs/${cond}
done

#autosort folders
for f in *;
do
	# baselines should always be C3H
    # will onnly use hets for Zic2
  if grep -q "wt" <<< "$f" && grep -q  "C3H"; then
		mv f baseline/inputs/baseline/
	elif grep -q "het" <<< "$f" && grep -q  "C3H"; then
		mv f mutants/inputs/mutants/
	elif grep -q "wt" <<< "$f" ; then
		mv f teatment/inputs/treatment/
	else
		mv f mut_treat/inputs/mut_treat/
	fi
done


















