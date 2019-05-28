"""
Make QC images of the registreres volumes
"""

from os.path import splitext, basename, join
from pathlib import Path

import SimpleITK as sitk
from logzero import logger as logging
import numpy as np

from lama import common
from lama.registration_pipeline.validate_config import LamaConfig


def make_qc_images_from_config(config: LamaConfig,
                               lama_specimen_dir: Path,
                               registerd_midslice_outdir: Path,
                               inverted_label_overlay_outdir: Path):
    """
    Generate midslice images for quick qc of registrations

    Parameters
    ----------
    config: The lama config
    lama_specimen_dir: The registration outdir. Should contain an 'output' folder
    registerd_midslice_outdir: Where to place the midslices
    inverted_label_overlay_outdir: Location of inverted labels

    Make qc images from:
        The final registration stage.
            What the volumes look like after registration
        The rigidly-registered images with the inverted labels overlaid
            This is a good indicator of regsitration accuracy

    """

    # Get the first registration stage (rigid)
    first_stage_id = config['registration_stage_params'][0]['stage_id']
    first_stage_dir = lama_specimen_dir / 'registrations' / first_stage_id

    # Get the final stage registration
    final_stage_id = config['registration_stage_params'][-1]['stage_id']
    final_stage_dir = lama_specimen_dir / 'registrations' / final_stage_id

    # Get the inverted labels dir, that will map onto the first stage registration
    inverted_label_id = config['registration_stage_params'][1]['stage_id']
    inverted_label_dir = lama_specimen_dir / 'inverted_labels' / inverted_label_id

    generate(first_stage_dir,
             final_stage_dir,
             inverted_label_dir,
             registerd_midslice_outdir,
             inverted_label_overlay_outdir)


def generate(first_stage_reg_dir: Path,
             final_stage_reg_dir: Path,
             inverted_labeldir: Path,
             out_dir_vols: Path,
             out_dir_labels: Path):
    """
    Generate a mid section slice to keep an eye on the registration process.

    When overalying the labelmaps, it depends on the registrated columes and inverted labelmaps being named identically
    """

    # Make midslice images of the final registered volumes
    if final_stage_reg_dir:
        p = common.get_file_paths(final_stage_reg_dir)
        if not p:
            logging.warn("Can't find output files for {}".format(final_stage_reg_dir))
            return

        for vol_path in p:

            # label_path = inverted_label_dir / vol_path.stem /vol_path.name

            vol_reader = common.LoadImage(vol_path)
            # label_reader = common.LoadImage(label_path)

            if not vol_reader:
                logging.error('error making qc image: {}'.format(vol_reader.error_msg))
                continue

            cast_img = sitk.Cast(sitk.RescaleIntensity(vol_reader.img), sitk.sitkUInt8)
            arr = sitk.GetArrayFromImage(cast_img)
            slice_ = np.flipud(arr[:, :, arr.shape[2] // 2])
            out_img = sitk.GetImageFromArray(slice_)

            base = splitext(basename(vol_reader.img_path))[0]
            out_path = join(out_dir_vols, base + '.png')
            sitk.WriteImage(out_img, out_path, True)

    # Make a sagittal mid-section overlay of inverted labels on input (rigidly-alingned) specimen
    if first_stage_reg_dir and inverted_labeldir:
        for vol_path in common.get_file_paths(first_stage_reg_dir):

            label_path = inverted_labeldir / vol_path.stem /vol_path.name
            vol_reader = common.LoadImage(vol_path)

            if not vol_reader:
                logging.error(f'cannnot create qc image from {vol_path}')
                return

            label_reader = common.LoadImage(label_path)
            if not label_reader:
                logging.error(f'cannot create qc image from label file {label_path}')
                return

            cast_img = sitk.Cast(sitk.RescaleIntensity(vol_reader.img), sitk.sitkUInt8)
            arr = sitk.GetArrayFromImage(cast_img)
            slice_ = np.flipud(arr[:, :, arr.shape[2] // 2])
            l_arr = label_reader.array
            l_slice_ = np.flipud(l_arr[:, :, l_arr.shape[2] // 2])

            base = splitext(basename(label_reader.img_path))[0]
            out_path = join(out_dir_labels, base + '.png')
            blend_8bit(slice_, l_slice_, out_path)


def blend_8bit(gray_img: np.ndarray, label_img: np.ndarray, out: Path, alpha: float=0.18):

    overlay_im = sitk.LabelOverlay(sitk.GetImageFromArray(gray_img),
                                   sitk.GetImageFromArray(label_img),
                                   alpha,
                                   0)
    sitk.WriteImage(overlay_im, out)


if __name__ == '__main__':
    config = Path('/mnt/IMPC_research/neil/E18.5/har_wt_test_080519/lama_e18_5.toml')

    root_dir = Path('/mnt/IMPC_research/neil/E18.5/har_wt_test_080519/output/baseline')
    inverted_label_dir = Path('/mnt/IMPC_research/neil/E18.5/har_wt_test_080519/QC/inverted_label_overlays')

    root_outdir = Path('/mnt/IMPC_research/neil/E18.5/har_wt_test_080519/QC')

    midslice_out = root_outdir / 'midslice'
    midslice_out.mkdir(exist_ok=True)
    inverted_label_dir = root_outdir / 'inverted_overlay'
    inverted_label_dir.mkdir(exist_ok=True)

    for spec_dir in root_dir.iterdir():
        spec_reg_dir  = spec_dir / 'output'

        c = LamaConfig(config)

        make_qc_images_from_config(c, spec_reg_dir, Path(midslice_out), Path(inverted_label_dir))

