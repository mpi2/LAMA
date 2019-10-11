from pathlib import Path
from lama.common import download_and_extract_zip


def main():
    url = 'https://www.mousephenotype.org/images/lama/test_data.zip'
    test_dir = Path(__file__).parent.parent / 'tests'

    download_and_extract_zip(url, test_dir)

