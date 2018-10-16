"""
Make QC images of the registreres volumes
"""

import SimpleITK as sitk
import common
from logzero import logger as logging
import numpy as np
from os.path import splitext, basename, join
from pathlib import Path
import matplotlib.pyplot as plt

def make_qc_images_from_config(config: dict, outdir: Path):
    """
    Parameters
    ----------
    config: The lama config used to get the corredt stage directories

    Make qc images from:
        The final registration stage.
            What the volumes look like after registration
        The rigidly-registered images with the inverted labels overlaid
            This is a good indicator of regsitration accuracy

    """

def make_qc_images(reg_stration_dir: Path,
                   inverted_label_dir: Path,
                   out_dir_vols: Path,
                   out_dir_labels):
    """
    Generate a mid section slice to keep an eye on the registration process.

    When overalying the labelmaps, it depends on the registrated columes and inverted labelmaps being named identically
    """
    for vol_path in common.get_file_paths(reg_stration_dir):

        label_path = inverted_label_dir / vol_path.stem /vol_path.name

        vol_reader = common.LoadImage(vol_path)
        label_reader = common.LoadImage(label_path)

        if not vol_reader:
            logging.error('error making qc image: {}'.format(vol_reader.error_msg))
            continue

        if not label_reader:
            logging.error('error making qc image. Cannot load label: {}'.format(label_reader.error_msg))

        cast_img = sitk.Cast(sitk.RescaleIntensity(vol_reader.img), sitk.sitkUInt8)
        arr = sitk.GetArrayFromImage(cast_img)
        slice_ = np.flipud(arr[:, :, arr.shape[2] // 2])
        out_img = sitk.GetImageFromArray(slice_)

        if label_reader:

            l_arr = label_reader.array
            l_slice_ = np.flipud(l_arr[:, :, l_arr.shape[2] // 2])

            base = splitext(basename(label_reader.img_path))[0]
            out_path = join(out_dir_labels, base + '.png')
            blend_8bit(slice_, l_slice_, out_path)

        # base = splitext(basename(vol_reader.img_path))[0]
        # out_path = join(out_dir_vols, base + '.png')
        # sitk.WriteImage(out_img, out_path, True)


def blend_8bit(gray_img: np.ndarray, label_img: np.ndarray, out: Path, alpha: float=0.18) -> np.ndarray:
    # blended = alpha * img1 + (1 - alpha) * img2
    overlay_im = sitk.LabelOverlay(sitk.GetImageFromArray(gray_img),
                                   sitk.GetImageFromArray(label_img),
                                   alpha,
                                   0)
    sitk.WriteImage(overlay_im, out)

if __name__ == '__main__':

    import shutil

    spec_root = Path('/home/neil/IMPC_research/neil/E14.5/mutants/specimens')
    for spec_dir in spec_root.iterdir():
        print(spec_dir.name)
        rigid_dir = spec_dir / 'output' / 'registrations' / 'rigid'
        label_dir = spec_dir / 'output' / 'inverted_labels' / 'similarity'

        qc_out_dir = spec_dir / 'output' / 'qc' / 'qc_inverted_labels_overlay'

        shutil.rmtree((qc_out_dir))
        qc_out_dir.mkdir(exist_ok=True)

        # in_lab_dir = Path('/home/neil/IMPC_research/neil/E14.5/mutants/specimens/1200014J11RIK/output/inverted_labels/similarity')
        # reg_dir = Path('/home/neil/IMPC_research/neil/E14.5/mutants/specimens/1200014J11RIK/output/registrations/rigid')
        # outlabls = Path('/home/neil/IMPC_research/neil/E14.5/mutants/specimens/1200014J11RIK/output/qc/test/labels')
        outvols = Path('/home/neil/IMPC_research/neil/E14.5/mutants/specimens/1200014J11RIK/output/qc/test/vols')

        make_qc_images(rigid_dir, label_dir, outvols , qc_out_dir)