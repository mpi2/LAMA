"""
Script to run the lama registration piepline
"""
import argparse

from lama.run_lama import RegistrationPipeline



if __name__ == "__main__":

    parser = argparse.ArgumentParser("The MRC Harwell image registration pipeline")
    parser.add_argument('-c', dest='config', help='Config file (YAML format)', required=True)
    args = parser.parse_args()

    RegistrationPipeline(args.config)
