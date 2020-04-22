# LAMA pipeline walkthroughs

See the lama wiki
The bash scripts in this folder show how to use LAMA for the following tasks.
The specimens used are from strain C57B/6N @ E14.5.

The scripts can be run from the command line or they can be used just as a guide.

First thing to do is download the data (~1.1GB) for running the walkthroughs. 
You can do this by running lama_get_walkthrough_data

1. Population average construction - population_average.sh
This script runs lama on 8 wild type males to create a population average.

2. Generation of data from wild type and mutant specimens - spatially_normalise_data.sh
The next step involves spatially normalising 8 wild type speciemns and 6 ACAN-/- mutant specimens. 
This is wehere the Jacobian determinants are also created along with the calulation of organ volumes.
This stage is not dependent on step 1 as it uses a premade population average and associated atlas.

3. Statistical analysis of the data from previous step - stats_with_BH_correction.sh
This step takes the spatially-normlised data from step 1. and does some statistical analysis on the organ volumes and Jacobian determinants.


