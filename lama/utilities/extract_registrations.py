
from pathlib import Path
import nrrd
import os
from lama import common

target_dir = Path("Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210722_Ku_conf_opt_v1/wild_type_and_mutant_data")

reg_paths = [spec_path for spec_path in common.get_file_paths(target_dir) if ('registrations' in str(spec_path))]

rigid_paths = [spec_path for spec_path in reg_paths if ('rigid' in str(spec_path))]

for path in rigid_paths:
    print(path)
    rigid, r_head = nrrd.read(path)

    print(str(os.path.basename(path)))
    nrrd.write(str(target_dir) + "/" + str(os.path.basename(path)), rigid)