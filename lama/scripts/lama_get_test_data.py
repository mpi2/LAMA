from pathlib import Path
from lama.common import download_and_extract_zip


def main():
    url = 'http://images.mousephenotype.org/lama/test_data.zip'
    test_dir = Path(__file__).parent.parent / 'tests'
    download_and_extract_zip(url, test_dir)


if __name__ == '__main__':
    main()