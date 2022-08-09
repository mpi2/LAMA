from pathlib import Path
import nrrd
import os
from lama import common
import SimpleITK as sitk
from scipy import ndimage
import numpy as np


def main(target_dir, labs_of_interest: list = [17]):
    label_paths = [spec_path for spec_path in common.get_file_paths(target_dir) if
                   ('inverted_labels' in str(spec_path))]

    rigid_paths = [rigid_path for rigid_path in common.get_file_paths(Path(target_dir / 'rigid'), extension_tuple=".nrrd")]


    rigid_paths.sort(key = lambda x: os.path.basename(x))
    label_paths.sort(key = lambda x: os.path.basename(x))

    #rigids = [nrrd.read(path) for path in rigid_paths]



    # Just get the label and write them
    for i, path in enumerate(label_paths):
        print(rigid_paths[i], path)

        label, l_head = nrrd.read(path)
        rigid, r_head = nrrd.read(rigid_paths[i])
        # print((~np.isin(label, labs_of_interest)))

        label[~np.isin(label, labs_of_interest)] = 0

        # get roi of label and rigid for scaling
        s = ndimage.find_objects(label)[-1]


        if isinstance(labs_of_interest, list):
            t = ndimage.find_objects(label)[int(min(labs_of_interest))]

            midpoint = [np.mean([t[0].start, s[0].stop]),
                        np.mean([t[1].start, s[1].stop]),
                        np.mean([t[2].start, s[2].stop])]
        else:
            midpoint = [np.mean([s[0].start, s[0].stop]),
                        np.mean([s[1].start, s[1].stop]),
                        np.mean([s[2].start, s[2].stop])]

        midpoint = [int(np.round(i)) for i in midpoint]

        p = 80


        crop_lab = label[midpoint[0] - p: midpoint[0] + p,
                   midpoint[1] - p: midpoint[1] + p,
                   midpoint[2] - p: midpoint[2] + p]



        crop_rig = rigid[midpoint[0] - p: midpoint[0] + p,
                   midpoint[1] - p: midpoint[1] + p,
                   midpoint[2] - p: midpoint[2] + p]

        # label[label not in labs_of_interest] = 0
        # print(str(os.path.basename(path)))
        os.makedirs(target_dir / "uncropped_labels", exist_ok=True)
        os.makedirs(target_dir / "cropped_labels", exist_ok=True)
        os.makedirs(target_dir / "cropped_rigids", exist_ok=True)

        file_name = target_dir / "uncropped_labels" / str(os.path.basename(path))

        cr_file_name = target_dir / "cropped_rigids" / str(os.path.basename(path))

        cl_file_name = target_dir / "cropped_labels" / str(os.path.basename(path))

        nrrd.write(str(file_name), label)
        nrrd.write(str(cr_file_name), crop_rig)
        nrrd.write(str(cl_file_name), crop_lab)
