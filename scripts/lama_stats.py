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
parser.add_argument('-c', '--config', dest='config', help='yaml config file contanign stats info', required=True)
args = parser.parse_args()
run(args.config)
