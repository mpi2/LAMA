"""
Upon completion of a a lama run, this scrpt will collate all the QC images into a single file format x to enable
the fast identification of issues

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


def run(root: Path, outdir,
        axial_slice=0,
        coronal_slice=0,
        sagittal_slice=0):
    """

    Parameters
    ----------
    root
        A Lama registrtion root for a line (or the baselines) containing multiple
    outpath
        Where to put the final image. Filetype is from extension (can use .png and .pdf at least)

    """
    # Create series of images specimen-based view

    d = DataIterator(root)

    spec: SpecimenDataPaths

    oris = [HtmlGrid('axial'),
            HtmlGrid('coronal'),
            HtmlGrid('sagittal')]

    for i, spec in enumerate(d):
        rc_qc_dir = spec.qc_red_cyan_dirs

        for grid in oris:
            spec.specimen_id
            spec_title = f'{spec.line_id}: {spec.specimen_id}'
            grid.next_row(title=spec_title)

            for img_path in natsorted((rc_qc_dir / grid.title).iterdir(), key=lambda x: x.stem): 
                relpath = Path(os.path.relpath(img_path, outdir))
                img_caption = f'{truncate_str(img_path.stem, 30)}'
                grid.next_image(relpath, img_caption)

    for grid in oris:
        ori_out = outdir / f'{grid.title}.html'
        grid.save(ori_out)


if __name__ =='__main__':
    import sys
    run(Path(sys.argv[1]), Path(sys.argv[2]))