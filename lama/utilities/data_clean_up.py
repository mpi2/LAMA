"""
LAMA produces lots of data. Sometimes we can can rid of much of it after wrads.
This script removes folde specified in a config file.

Work in progress


Current example toml
--------------------
folders_to_rm = [
'resolution_images']
--------------------

This will search all suddirectories an delete any folder called 'resolution_images'
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

    for dir_ in config['folders_to_rm']:
        rm_by_name(root_dir, dir_)


if __name__ == '__main__':
    import sys
    config_path = sys.argv[1]
    root_dir = sys.argv[2]
    run(config_path, Path(root_dir))
