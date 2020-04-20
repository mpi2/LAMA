"""
Download and unzip the data for the walk-throughs
"""

from pathlib import Path
import sys


from lama.common import download_and_extract_zip


def main():
    url = 'http://images.mousephenotype.org/lama/walkthrough_data.zip'


    if len(sys.argv) > 1:
        unzip_dir = sys.argv[1]
    else:
        unzip_dir = Path(__file__).parent

    download_and_extract_zip(url, unzip_dir)


if __name__ == '__main__':
    main()
