#!/bin/bash

# Complete image processing script created by Amrit and Kyle


# cropping scans
mv cropper.py nrrd_out/

ls nrrd_out

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





















