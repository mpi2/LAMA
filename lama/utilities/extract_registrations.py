
from pathlib import Path
import nrrd
import os
from lama import common




def main(target_dir):


    reg_paths = [spec_path for spec_path in common.get_file_paths(target_dir) if ('registrations' in str(spec_path))]

    rigid_paths = [spec_path for spec_path in reg_paths if ('rigid' in str(spec_path))]




    rigid_paths.sort()

    for path in rigid_paths:
        rigid, r_head = nrrd.read(path)
        os.makedirs(target_dir / "rigid", exist_ok=True)

        file_name = target_dir / "rigid" / str(os.path.basename(path))
        nrrd.write(str(file_name), rigid)