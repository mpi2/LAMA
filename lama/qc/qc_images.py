"""
Make QC images of the registered volumes
"""

from os.path import splitext, basename, join
from pathlib import Path

import SimpleITK as sitk
from logzero import logger as logging
import numpy as np
from skimage.exposure import rescale_intensity
from skimage.io import imsave

from lama import common
from lama.registration_pipeline.validate_config import LamaConfig
from lama.elastix import IGNORE_FOLDER

MOVING_RANGE = (0, 180)  # Rescale the moving image to these values for the cyan/red overlay


def make_qc_images_from_config(config: LamaConfig,
                               lama_specimen_dir: Path,
                               registerd_midslice_outdir: Path,
                               inverted_label_overlay_outdir: Path,
                               cyan_red_dir: Path,
                               target: Path):
    """
    Generate mid-slice images for quick qc of registrations

    Parameters
    ----------
    config
        The lama config
    lama_specimen_dir
        The registration outdir. Should contain an 'output' folder
    registerd_midslice_outdir
        Where to place the midslices
    inverted_label_overlay_outdir
        Location of inverted labels
    cyan_red_dir
        Where to put target(cyan) moving(red) overlays
    target
        The tareget image to dispaly in cyan

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
             inverted_label_overlay_outdir,
             cyan_red_dir,
             target)


def generate(first_stage_reg_dir: Path,
             final_stage_reg_dir: Path,
             inverted_labeldir: Path,
             out_dir_vols: Path,
             out_dir_labels: Path,
             out_dir_cyan_red: Path,
             target: Path):
    """
    Generate a mid section slice to keep an eye on the registration process.

    When overalying the labelmaps, it depends on the registrated columes and inverted labelmaps being named identically
    """

    target = common.LoadImage(target).array
    target_slice = target[:, :, target.shape[2] // 2]

    # Make mid-slice images of the final registered volumes
    if final_stage_reg_dir:
        p = common.get_file_paths(final_stage_reg_dir, ignore_folder=IGNORE_FOLDER)
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

            rc_slice = overlay_cyan_red(target_slice, slice_)
            imsave(out_dir_cyan_red / (base + '.png'), rc_slice)

    # Make a sagittal mid-section overlay of inverted labels on input (rigidly-aligned) specimen
    if first_stage_reg_dir and inverted_labeldir:
        for vol_path in common.get_file_paths(first_stage_reg_dir, ignore_folder=IGNORE_FOLDER):

            vol_reader = common.LoadImage(vol_path)

            if not vol_reader:
                logging.error(f'cannnot create qc image from {vol_path}')
                return

            label_path = inverted_labeldir / vol_path.stem / vol_path.name

            if label_path.is_file():
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
            else:
                logging.info('No inverted label found. Skipping creation of inverted label-image overlay')


def blend_8bit(gray_img: np.ndarray, label_img: np.ndarray, out: Path, alpha: float=0.18):

    overlay_im = sitk.LabelOverlay(sitk.GetImageFromArray(gray_img),
                                   sitk.GetImageFromArray(label_img),
                                   alpha,
                                   0)
    sitk.WriteImage(overlay_im, out)


def overlay_cyan_red(target: np.ndarray, specimen: np.ndarray) -> np.ndarray:
    """
    Create a cyan red overlay
    Parameters
    ----------
    target
    specimen`

    Returns
    -------

    """
    if target.shape != specimen.shape:
        raise ValueError('target and specimen must be same shape')

    specimen = np.clip(specimen, 0, 255)
    specimen = rescale_intensity(specimen, out_range=(MOVING_RANGE))

    rgb = np.zeros([*target.shape, 3], np.uint8)

    rgb[..., 0] = target
    rgb[..., 1] = specimen
    rgb[..., 2] = specimen

    return np.flipud(rgb)


