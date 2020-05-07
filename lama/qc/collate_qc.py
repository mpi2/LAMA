"""
Upon completion of a a lama run, this scrpt will collate all the QC images into a single file format x to enable
the fast identification of issues

"""
from pathlib import Path
from math import ceil

from skimage import io
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

from lama.paths import DataIterator, SpecimenDataPaths


def run(root: Path, outpath, num_cols: int=6):
    """

    Parameters
    ----------
    root
        A Lama registrtion root for a line (or the nbaselines) containing multiple
    outpath
        Where to put the final image. Filetype is from extension (can use .png and .pdf at least)
    num_cols

    """
    d = DataIterator(root)

    num_specimens = len(d)
    num_rows = ceil(num_specimens / num_cols)

    fig_height = 100 * num_rows
    fig = plt.figure(figsize=(20., fig_height))

    grid = ImageGrid(fig, 111,  # similar to subplot(111)
                     nrows_ncols=(num_rows, num_cols),
                     axes_pad=[0.1, 0.5],  # pad between axes in inch.
                     )
    s: SpecimenDataPaths

    for ax, s in zip(grid, d):
        # Lets get the cyan red overlays
        img = io.imread(s.qc_red_cyan)
        ax.imshow(img)
        ax.set_title(f'{s.line_id}\n{s.specimen_id}')
    fig.savefig(outpath, bbox_inches='tight')




if __name__ =='__main__':
    run(Path('/mnt/IMPC_media/LAMA_staging/e15_5/240420_initial_test/centres/output_bu'), '/home/neil/Desktop/t/grid.png')