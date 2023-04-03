from pathlib import Path
import nrrd
import SimpleITK as sitk
from scipy import ndimage



img_path = '/mnt/bit_nfs/neil/impc_e15_5/phenotyping_tests/BCM/stability_101120/output/baseline/output/baseline/2086980_download/output/registrations/rigid/2086980_download/2086980_download.nrrd'
labels_path = '/mnt/bit_nfs/neil/impc_e15_5/phenotyping_tests/BCM/stability_101120/output/baseline/output/baseline/2086980_download/output/inverted_labels/similarity/2086980_download/2086980_download.nrrd'
outdir = Path('/home/neil/Desktop/t/cluster')

label_of_interest = 10 # medial liver lobe

img, i_head = nrrd.read(img_path)
labels, _ = nrrd.read(labels_path)

# Get the bounding box of the inverted label
labels[labels != label_of_interest] = 0

s = ndimage.find_objects(labels)[9]

# Add some padding
p = 30

# Get a rough ROI based on the lama segmentation
roi = img[s[0].start - p: s[0].stop + p ,
       s[1].start -p: s[1].stop + p,
       s[2].start -p: s[2].stop + p]

outpath = outdir / 'test_roi.nrrd'
nrrd.write(str(outpath), roi, header=i_head)
