from pathlib import Path
import urllib, io
import urllib.request
import zipfile

url = 'https://www.mousephenotype.org/images/lama/test_data.zip'

current_dir = Path(__file__).parent

test_data_path = current_dir.resolve()

print('Downloading test data')

remotezip = urllib.request.urlopen(url)
zipinmemory = io.BytesIO(remotezip.read())
zip = zipfile.ZipFile(zipinmemory)
zip.extractall(test_data_path)

print(f'Test data downloaded and extracted to {test_data_path}')
