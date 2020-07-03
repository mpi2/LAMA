"""
Remove specimen output folders out of the main output directory.
This might be needed if, for example, QC issues were identified in a specimen which then needed removing before further
analysis
"""
from lama.paths import specimen_iterator
from typing import List
from lama.common import csv_read_lines
from pathlib import Path
import shutil


def move_stuff(root_dir: Path, out_dir: Path, spec_ids: [List]):
    for line_dir, spec_dir in specimen_iterator(root_dir):
        if spec_dir.name in spec_ids:
            print(f"Moving {spec_dir.name}")

            dest = out_dir / spec_dir.name
            shutil.move(spec_dir, dest)


if __name__ == '__main__':
    import sys
    _root_dir = Path(sys.argv[1])
    _out_dir = Path(sys.argv[2])
    _ids_file = Path(sys.argv[3])

    _ids = csv_read_lines(_ids_file)

    move_stuff(_root_dir, _out_dir, _ids)