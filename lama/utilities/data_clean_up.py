"""
LAMA produces lots of data. Sometimes we can get rid of much of it afterwrads.
This script removes folders specified in a config file.

This is a work in progress


Current example toml
--------------------
folders_to_rm = [
'resolution_images',
'QC']
--------------------

python3 data_clean_up.py config root_dir

This will recursively search directories and delete any folder called in the list
"""

from pathlib import Path
import shutil
import toml


def rm_by_name(root: Path, name: str):
    dirs = root.glob(f'**/{name}')
    for d in dirs:
        shutil.rmtree(d)


def run(config_path: str, root_dir: Path):

    config = toml.load(config_path)
    print(f"deleting {config['folders_to_rm']}")

    for dir_ in config['folders_to_rm']:
        rm_by_name(root_dir, dir_)


if __name__ == '__main__':
    import sys
    config_path = sys.argv[1]
    root_dir = sys.argv[2]
    run(config_path, Path(root_dir))
