"""
Upon completion of a a lama run, this scrpt will collate all the QC images into a single html to enable
the rapid identification of issues.

"""
from pathlib import Path
import os
from math import ceil
from matplotlib import figure
from skimage import io
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from collections import defaultdict
import addict
from lama.common import truncate_str
from lama.qc.img_grid import HtmlGrid
from natsort import natsorted

from lama.paths import DataIterator, SpecimenDataPaths


def make_grid(root: Path, outdir, qc_type='red_cyan', height='auto'):
    """

    Parameters
    ----------
    root
        A Lama registrtion root for a line (or the baselines) containing multiple specimens
    outpath
        Where to put the final image. Filetype is from extension (can use .png and .pdf at least)
    height:
        the css height property for the img (Use 'auto' or px. Percentage sclaing messes things up

    """
    # Create series of images specimen-based view

    d = DataIterator(root)

    spec: SpecimenDataPaths

    oris = [HtmlGrid('axial'),
            HtmlGrid('coronal'),
            HtmlGrid('sagittal')]

    for i, spec in enumerate(d):

        try:
            spec.setup()
        except FileNotFoundError as e:
            print(f'Skipping {spec.specimen_id}\n{e}')
            continue

        if qc_type == 'red_cyan':
            rc_qc_dir = spec.qc_red_cyan_dirs
        elif qc_type == 'grey':
            rc_qc_dir = spec.qc_grey_dirs
        elif qc_type == 'labels':
            rc_qc_dir = spec.qc_inverted_labels

        for grid in oris:
            spec.specimen_id
            spec_title = f'{spec.line_id}: {spec.specimen_id}'
            grid.next_row(title=spec_title)


            # s = list((rc_qc_dir / grid.title).iterdir())
            for img_path in natsorted((rc_qc_dir / grid.title).iterdir(), key=lambda x: x.stem):
                relpath = Path(os.path.relpath(img_path, outdir))
                img_caption = f'{truncate_str(img_path.stem, 30)}'
                tooltip = f'{spec.line_id}:{spec.specimen_id}:{img_path.stem}'
                grid.next_image(relpath, img_caption, tooltip)

        for grid in oris:
            ori_out = outdir / f'{grid.title}.html'
            grid.save(ori_out, height)


def run(reg_root: Path, out_root: Path, height):

    rc_dir = out_root / 'red_cyan'
    rc_dir.mkdir(exist_ok=True)
    make_grid(reg_root,  rc_dir, 'red_cyan', height=height)
    #
    g_dir = out_root / 'greyscales'
    g_dir.mkdir(exist_ok=True)
    make_grid(reg_root, g_dir, 'grey', height=height)

    g_dir = out_root / 'inverted_labels'
    g_dir.mkdir(exist_ok=True)

    try:
        make_grid(reg_root, g_dir, 'labels', height)
    except FileNotFoundError:
        print('Cannot find inverted label overlays. Skipping')


if __name__ =='__main__':

    import argparse
    parser = argparse.ArgumentParser("Genate HTML registration image reports")

    parser.add_argument('-i', '--indir', dest='indir', help='A lama registration output directory containing one or more line directories',
                        required=True)
    parser.add_argument('-o', '--out', dest='out', help='output directory',
                        required=True)
    parser.add_argument('-height', '--height', dest='height', help='The height of the images. eg "auto", 200px',
                        required=False, default='auto')

    args = parser.parse_args()
    run(Path(args.indir), Path(args.out), args.height)