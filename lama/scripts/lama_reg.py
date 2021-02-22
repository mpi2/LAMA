"""
Script to run the lama registration piepline.

Note: This does not work proerly any more.
"""
import argparse
from pathlib import Path

from lama.registration_pipeline import run_lama


def main():
    """
    setup tools console_scripts used to allow command lien running of scripts needs to have a function that takes
    zero argumnets.
    """
    from lama.elastix import invert_transforms
    parser = argparse.ArgumentParser("The LAMA registration pipeline")
    parser.add_argument('-c', dest='config', help='Config file (TOML format)', required=True)
    args = parser.parse_args()

    run_lama.run(Path(args.config))


if __name__ == "__main__":
    main()
