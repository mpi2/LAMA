"""
Iterate over qc images from a lama run
For each specimen
    show a series of qc images
    in terminal asccept options such as f(fail), p(pass), c(check)
    After option then add optional comment
    Add row to a datframe
    Save dataframe

Example
-------
collate_qc_new.py output/baseline/output qc_status.csv

"""
from pathlib import Path
from skimage.io import imread
from natsort import natsorted
from itertools import chain
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from screeninfo import get_monitors
from mpl_toolkits.axes_grid1 import ImageGrid
from skimage.transform import resize

from lama.paths import DataIterator, SpecimenDataPaths


def show_(fig, img_path, row, ncol, idx, cmap=None, flip=False):
    fig.add_subplot(row, ncol, idx)
    img = imread(img_path)
    if flip:
        img = np.flipud(img)
    plt.imshow(img, cmap=cmap)
    plt.gca().set_axis_off()
    # plt.subplots_adjust(top=1, bottom=0, right=1, left=0,
    #                     hspace=0, wspace=0)
    plt.margins(0, 0)


def do_qc(root: Path, outpath: Path):
    """

    Parameters
    ----------
    root
        A Lama registrtion root for a line (or the baselines) containing multiple specimens
    outpath
        Qc csv

    """
    # Create series of images specimen-based view
    # mpl.rcParams['figure.raise_window']
    d = DataIterator(root)

    spec: SpecimenDataPaths

    status_map = {
        'f': 'failed',
        'p': 'passed',
        'c': 'check'
    }

    plt.ion()



    if outpath.is_file():
        df = pd.read_csv(outpath, index_col=0)
    else:
        df = pd.DataFrame(columns=['line', 'spec_id', 'status', 'comment'])

    # get first monitor. Could also look for largest monitor?
    mon = get_monitors()[0]
    width = mon.width
    width = width - (width * 0.25)  # Save som sapce for the console
    height = mon.height
    dpi = 100
    # What size does the figure need to be in inches to fit the image?
    figsize = width / float(dpi), height / float(dpi)
    f = plt.figure(figsize=figsize, dpi=dpi)

    grid = ImageGrid(f, 111,  # similar to subplot(111)
                     nrows_ncols=(2, 8),  # creates 2x2 grid of axes
                     axes_pad=0.1,  # pad between axes in inch.
                     )

    for i, spec in enumerate(d):

        try:
            spec.setup()
        except FileNotFoundError as e:
            print(f'Skipping {spec.specimen_id}\n{e}')
            continue

        if spec.specimen_id in df.spec_id.values:
            print(f'skipping {spec.specimen_id}')
            continue

        # Get the overlays. Use natsort to sort asceding by the slice number file suffix
        ol_dir = spec.qc_inverted_labels
        ax_ol = natsorted([x for x in (ol_dir / 'axial').iterdir()])
        sa_ol = natsorted([x for x in (ol_dir / 'sagittal').iterdir()])
        co_ol = natsorted([x for x in (ol_dir / 'coronal').iterdir()])

        reg_dir = spec.qc_grey_dirs
        ax_ol_r = natsorted((reg_dir / 'axial').iterdir(), key= lambda x: x.name)[-1]
        sa_ol_r = natsorted((reg_dir / 'sagittal').iterdir() ,key= lambda x: x.name)[-1]
        co_ol_r = natsorted((reg_dir / 'coronal').iterdir(), key= lambda x: x.name)[-1]


        all_p = chain(ax_ol, sa_ol, co_ol,  [ax_ol_r, sa_ol_r, co_ol_r])
        # num_cols = len(ax_ol) *2

        len_p = len(list(chain(ax_ol, sa_ol, co_ol)))
        i = 0
        for ax, im in zip(grid, all_p ):
            i+=1
            img = imread(im)
            # new_shape = [x*2for x in img.shape[:2]]
            # new_shape.append(3)
            # img = resize(img, new_shape, order=2)
            if i < len_p:
                ax.imshow(img)
            else:
                ax.imshow(img, cmap='gray')

        f.suptitle(f'{spec.line_id} - {spec.specimen_id}')

        plt.show()
        plt.pause(0.1)
        while True:
            status = input("Enter status: ")
            if status not in status_map.keys():
                print(f'status must be one of {status_map.keys()}')
            else:
                break

        comment = input('Optional comment:')
        df.loc[len(df)] = [spec.line_id, spec.specimen_id, status, comment]
        df.to_csv(outpath)
        plt.close()



if __name__ =='__main__':
    import sys
    do_qc(Path(sys.argv[1]), Path(sys.argv[2]))