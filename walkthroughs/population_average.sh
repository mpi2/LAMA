#!/bin/bash

#### Some preliminary stuff to get it working ####
# Get the current working directory of this script
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# printf this beforwe each command to get indentation and dollar prompt
c="\n\t$ "

# Print out the commands except echos and reads
trap '! [[ "$BASH_COMMAND" =~ ^(echo|read|printf) ]] && \
cmd=`eval echo "$BASH_COMMAND" 2>/dev/null` && echo "$cmd"' DEBUG



#### Start the walk-through ####
printf $c
cd "$parent_path"/data/population_average

# Check to see if all images are same size and bit depth and then pad 
read -p "Press enter to get image metadata:"
printf $c
lama_img_info -i inputs

echo "As we can see, the dimensions (z, y, x) are not all the same."
read -p "Press enter to pad all volumes in place:"

# Pad volumes inplace
printf $c
lama_pad_volumes -i inputs/ --clobber

# Check the sizes agains
lama_img_info -i inputs
echo "Now all the volumes should be the same size:"

# Run the population average construction
read -p "Press enter to run the registration:"
printf $c
lama_reg -c pop_avg.toml

echo -e "Registration finished!\n"

# The resultant population average should now be in tuutorials/population_average/output/Averages
# Now see *wild_type_and_mutant_data*


