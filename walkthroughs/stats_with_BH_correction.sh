#!/bin/bash

# We've now speatially normalized our wild tpye specimen controls as well as specimens from a single knockout line.
# Now we can run the statistical analysis on the organ_volumes and Jacobian determinants.

# Multiple testing correction
# As we have only one knockout line, we will perform Benjamini-Hochburg FDR correction across the organ volume or
# Jacobian determinant p-values. 


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

# We will have one output directory per knockout line. Which is one in this case.
read -p "Press enter to run the stats pipeline"
lama_stats -c stats.toml -w ../wild_type_and_mutant_data/baseline  -m ../wild_type_and_mutant_data/mutants \
-o stats_output -t ../wild_type_and_mutant_data/target



