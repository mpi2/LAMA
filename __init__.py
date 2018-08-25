import logging
import subprocess
import os
import warnings


print("LAMA phenotype detection package\n\n")

# When running from a docker image, get much warnings from sklearn, pandas etc. Turn it off
warnings.filterwarnings("ignore")
