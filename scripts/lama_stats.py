#! /usr/bin/env python3

"""
Usage
-----
This script should be run from the root directory (one up from this script)

$~ scripts/lama_stats.py -c <path to stats config>
"""


import sys
import argparse
from lama import common
from lama.stats.run_lama_stats import run

sys.excepthook = common.excepthook_overide

parser = argparse.ArgumentParser("Stats component of the phenotype detection pipeline")
parser.add_argument('-w', '--wildtype_dir', dest='wt_dir', help='wild egistration output root directory', required=True)
parser.add_argument('-m', '--mutant_dir', dest='mut_dir', help='mutant registration output root directory', required=True)
parser.add_argument('-o', '--output_dir', dest='out_dir', help='Output directory to put collated results from all lines', required=True)
parser.add_argument('-m', '--mutant_dir', dest='mut_dir', help='mutant registration output root directory', required=True)
args = parser.parse_args()
run(args.config)
