"""
Make QC images of the registered volumes
"""

from os.path import splitext, basename, join
from pathlib import Path
from typing import List

import SimpleITK as sitk
from logzero import logger as logging
import numpy as np
from skimage.exposure import rescale_intensity, match_histograms
from skimage.io import imsave
from skimage.measure import regionprops

from lama import common
from lama.elastix import RESOLUTION_IMGS_DIR, IMG_PYRAMID_DIR
from lama.paths import LamaSpecimenData

INTENSITY_RANGE = (0, 255)  # Rescale the moving image to these values for the cyan/red overlay


def make_qc_images(lama_specimen_dir: Path,
                   target: Path,
                   outdir: Path,
                   mask: Path,
                   reverse_reg_propagation: bool = False):
    """
    Generate mid-slice images for quick qc of registration process.

    Parameters
    ----------
    lama_specimen_dir
        The registration outdir. Should contain an 'output' folder
    target
        The target image to display in cyan
    outdir
        Where to put the qc images
    mask
        Used to identify the embryo in the image so we can display useful info
    reverse_reg_propagation
        Whether to overlay on orginal unregistered input (False) or on initial, probably rigid, registered image (True)

    Notes
    -----
    Make qc images from:
        The final registration stage.
            What the volumes look like after registration
        The rigidly-registered images with the inverted labels overlaid
            This is a good indicator of regsitration accuracy

    """
    target = common.LoadImage(target).array
    # Make qc images for all stages of registration including any resolution images
    try:
        paths = LamaSpecimenData(lama_specimen_dir).setup()
    except FileNotFoundError as e:
        logging.exception(f'cannot find specimen directory\n{e}')

    # Order output dirs by qc type
    red_cyan_dir = outdir / 'red_cyan_overlays'
    greyscale_dir = outdir / 'greyscales'
    red_cyan_dir.mkdir(exist_ok=True)
    greyscale_dir.mkdir(exist_ok=True)

    try:
        for i, (stage, img_path) in enumerate(paths.registration_imgs()):
            img = common.LoadImage(img_path).array
            _make_red_cyan_qc_images(target, img, red_cyan_dir, greyscale_dir, img_path.stem, i, stage)

        if paths.inverted_labels_dirs:

            # First reg img will the rigid-registered image
            first_reg_dir = paths.reg_dirs[0]

            if reverse_reg_propagation:
                # We have a reverse registration method of label propagation so we overlay the labels that were transformed
                # using the reverse registrtion transform (the final defoemable stage) as the target will have been the
                # Rigid input
                inverted_label_dir = paths.inverted_labels_dirs[-1]
            else:
                # The labels were propagated using the inverse transfrom method. Therefore we overlay the labels transformed
                # using the tforms up to the inverted affine stage onto the rigid input.
                # (could do inverted rigid labels overalid on orginal input, but on rigid allllows us to compare specimens
                # more easily using this method)
                inverted_label_dir = paths.inverted_labels_dirs[-2]  # -2 should be the step after rigid

            inverted_label_overlays_dir = outdir / 'inverted_label_overlay'
            inverted_label_overlays_dir.mkdir(exist_ok=True)

            _overlay_labels(first_reg_dir,
                            inverted_label_dir,
                            inverted_label_overlays_dir,
                            mask=mask)
    except FileNotFoundError: # 220221 bodge. lama_reg creates a different file structure tha job_runner. Need to harmonise
        logging.error('No QC images made. This maybe because you used lama_reg rather than lama_job_runner')


def _overlay_labels(first_stage_reg_dir: Path,
                    inverted_labeldir: Path,
                    out_dir_labels: Path,
                    mask: Path=None):
    """
    Overlay the first registrated image (rigid) with the corresponding inverted labels
    It depends on the registered volumes and inverted label maps being named identically
    """
    if mask:
        mask = sitk.GetArrayFromImage(sitk.ReadImage(str(mask)))
        rp = regionprops(mask)
        # Get the largest label. Likley only one from the mask
        mask_props = list(reversed(sorted(rp, key=lambda x: x.area)))[0]
        bbox = mask_props['bbox']

    for vol_path in common.get_file_paths(first_stage_reg_dir, ignore_folders=[RESOLUTION_IMGS_DIR, IMG_PYRAMID_DIR]):

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
            base = splitext(basename(label_reader.img_path))[0]
            l_arr = label_reader.array

            def sag(idx_):
                slice_sag = np.flipud(arr[:, :, idx_])
                l_slice_sag = np.flipud(l_arr[:, :, idx_])
                sag_dir = out_dir_labels / 'sagittal'
                sag_dir.mkdir(exist_ok=True)
                out_path_sag = sag_dir / f'{base}_{idx_}.png'
                _blend_8bit(slice_sag, l_slice_sag, out_path_sag)
            if mask is None: # get a few slices from middle
                sag_indxs = np.linspace(0, arr.shape[2], 8, dtype=np.int)[2:-2]
            else:
                sag_start = bbox[2]
                sag_end = bbox[5]
                sag_indxs = np.linspace(sag_start, sag_end, 6, dtype=np.int)[1:-1] # Take the 4 inner slices
            for idx in sag_indxs:
                sag(idx)

            def ax(idx_):
                slice_ax = arr[idx_, :, :]
                l_slice_ax = l_arr[idx_, :, :]
                ax_dir = out_dir_labels / 'axial'
                ax_dir.mkdir(exist_ok=True)
                out_path_ax = ax_dir / f'{base}_{idx_}.png'
                _blend_8bit(slice_ax, l_slice_ax, out_path_ax)
            if mask is None: # get a few slices from middle
                ax_indxs = np.linspace(0, arr.shape[0], 8, dtype=np.int)[2:-2]
            else:
                ax_start = bbox[0]
                ax_end = bbox[3]
                ax_indxs = np.linspace(ax_start, ax_end, 6, dtype=np.int)[1:-1]
            for idx in ax_indxs:
                ax(idx)

            def cor(idx_):
                slice_cor = np.flipud(arr[:, idx_, :])
                l_slice_cor = np.flipud(l_arr[:, idx_, :])
                cor_dir = out_dir_labels / 'coronal'
                cor_dir.mkdir(exist_ok=True)
                out_path_cor = cor_dir / f'{base}_{idx_}.png'
                _blend_8bit(slice_cor, l_slice_cor, out_path_cor)
            if mask is None: # get a few slices from middle
                cor_indxs = np.linspace(0, arr.shape[1], 8, dtype=np.int)[2:-2]
            else:
                cor_start = bbox[1]
                cor_end = bbox[4]
                cor_indxs = np.linspace(cor_start, cor_end, 6, dtype=np.int)[1:-1]
            for idx in cor_indxs:
                cor(idx)

        else:
            logging.info('No inverted label found. Skipping creation of inverted label-image overlay')


def _blend_8bit(gray_img: np.ndarray, label_img: np.ndarray, out: Path, alpha: float=0.18):

    overlay_im = sitk.LabelOverlay(sitk.GetImageFromArray(gray_img),
                                   sitk.GetImageFromArray(label_img),
                                   alpha,
                                   0)
    sitk.WriteImage(overlay_im, str(out))


def _red_cyan_overlay(slice_1, slice_2) -> np.ndarray:

    rgb = np.zeros([*slice_1.shape, 3], np.uint8)
    rgb[..., 0] = slice_1
    rgb[..., 1] = slice_2
    rgb[..., 2] = slice_2
    return rgb


def _make_red_cyan_qc_images(target: np.ndarray,
                             specimen: np.ndarray,
                             out_dir: Path,
                             grey_cale_dir: Path,
                             name: str,
                             img_num: int,
                             stage_id: str) -> List[np.ndarray]:
    """
    Create a cyan red overlay
    Parameters
    ----------
    target
    specimen`
    img_num
        A number to prefix onto the qc image so that when browing a folder the images will be sorteed
    Returns
    -------
    0: axial
    1: coronal
    2: sagittal
    """
    oris = ['axial', 'coronal', 'sagittal']

    def get_ori_dirs(root: Path):
        res = []
        for ori_name in oris:
            dir_ = root / ori_name
            dir_.mkdir(exist_ok=True)
            res.append([dir_, ori_name])
        return res

    if target.shape != specimen.shape:
        raise ValueError('target and specimen must be same shape')

    # specimen = np.clip(specimen, 0, 255)
    specimen = rescale_intensity(specimen, out_range=INTENSITY_RANGE).astype(np.uint8)
    target = rescale_intensity(target, out_range=INTENSITY_RANGE).astype(np.uint8)

    def get_slices(img):
        slices = []
        for ax in [0,1,2]:
            slices.append(np.take(img, img.shape[ax] // 2, axis=ax))
        return slices

    t = get_slices(target)
    s = get_slices(specimen)

    # histogram match the specimen to the target for each orientation slice
    #   This produces bad results sometimes. Move to adaptive histogram equalization

    # Todo: histogram matchng swtitched off as it makes the specimen disapaer. Fix this
    # s = [exposure.equalize_adapthist(x, clip_limit=0.03) for x in s]
    # Try mean normalisation
    # for ti, si in zip(t, s):
    #     med_t = ti[ti > 5].mean()
    #     med_s = si[si > 5].mean()
    #     diff = med_t - med_s
    #     si += int(diff)
    s = [match_histograms(si, ti) for si, ti in zip(s, t)]

    rc_oris = get_ori_dirs(out_dir)
    grey_oris = get_ori_dirs(grey_cale_dir)

    # put slices in folders by orientation
    for i in range(len(oris)):
        grey = s[i]
        rgb = _red_cyan_overlay(s[i], t[i])
        if not oris[i] == 'axial':
            rgb = np.flipud(rgb)
            grey = np.flipud(grey)
        imsave(rc_oris[i][0] / f'{img_num}_{stage_id}_{name}_{rc_oris[i][1]}.png', rgb)
        imsave(grey_oris[i][0]  / f'{img_num}_{stage_id}_{name}_{rc_oris[i][1]}.png', grey)

