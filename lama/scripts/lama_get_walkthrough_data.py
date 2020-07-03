"""
Download and unzip the data for the walk-throughs
"""

from pathlib import Path
import sys


from lama.common import download_and_extract_zip


def main():
    url = 'http://images.mousephenotype.org/lama/lama_walkthroughs.zip'

    # Unzip into cwd
    unzip_dir = Path().cwd()

    download_and_extract_zip(url, unzip_dir)


if __name__ == '__main__':
    main()
