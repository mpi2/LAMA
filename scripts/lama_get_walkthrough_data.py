"""
Download and unzip the data for the walk-throughs
"""

from pathlib import Path
import sys


from lama.common import download_and_extract_zip


def main():
    url = 'http://images.mousephenotype.org/lama/tutorial_data.zip'

    if len(sys.argv) < 2:
        print(f'Usage: {Path(sys.argv[0]).name} directory/to/extract_walkthrough_data_to')
        return

    unzip_dir = Path(sys.argv[1])
    if not unzip_dir.is_dir():
        print(f"{unzip_dir} is not a directory")
        return

    download_and_extract_zip(url, unzip_dir)


if __name__ == '__main__':
    main()