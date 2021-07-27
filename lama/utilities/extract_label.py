from pathlib import Path
import nrrd
import os
from lama import common
import SimpleITK as sitk
from scipy import ndimage

# img_path = '/mnt/bit_nfs/neil/impc_e15_5/phenotyping_tests/BCM/stability_101120/output/baseline/output/baseline/2086980_download/output/registrations/rigid/2086980_download/2086980_download.nrrd'
# labels_path = '/mnt/bit_nfs/neil/impc_e15_5/phenotyping_tests/BCM/stability_101120/output/baseline/output/baseline/2086980_download/output/inverted_labels/similarity/2086980_download/2086980_download.nrrd'
# outdir = Path('/home/neil/Desktop/t/cluster')

label_of_interest = 17  # lateral brain ventricle

target_dir = Path("Z:/ArkellLab/Lab Members/Kyle/PhD/vmshare/Zic2_Kumba_LAMA/210722_Ku_conf_opt_v1/wild_type_and_mutant_data")

label_paths = [spec_path for spec_path in common.get_file_paths(target_dir) if ('inverted_labels' in str(spec_path))]


# Just get the label and write them
for path in label_paths:
    print(path)
    label, l_head = nrrd.read(path)
    label[label != label_of_interest] = 0
    print(str(os.path.basename(path)))
    nrrd.write(str(target_dir) + "/" + str(os.path.basename(path)), label)