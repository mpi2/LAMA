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
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from lama.paths import DataIterator, SpecimenDataPaths


def show_(fig, img_path, ncol, idx, cmap=None, flip=False):
    fig.add_subplot(1, ncol, idx)
    img = imread(img_path)
    if flip:
        img = np.flipud(img)
    plt.imshow(img, cmap=cmap)
    plt.gca().set_axis_off()
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0,
                        hspace=0, wspace=0)
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

    d = DataIterator(root)

    spec: SpecimenDataPaths

    status_map = {
        'f': 'failed',
        'p': 'passed',
        'c': 'check'
    }

    n_cols = 6

    plt.ion()

    if outpath.is_file():
        df = pd.read_csv(outpath, index_col=0)
    else:
        df = pd.DataFrame(columns=['line', 'spec_id', 'status', 'comment'])

    for i, spec in enumerate(d):

        try:
            spec.setup()
        except FileNotFoundError as e:
            print(f'Skipping {spec.specimen_id}\n{e}')
            continue

        if spec.specimen_id in df.spec_id.values:
            print(f'skipping {spec.specimen_id}')
            continue

        ol_dir = spec.qc_inverted_labels
        ax_ol = ol_dir / 'axial' / f'{spec.specimen_id}.png'
        sa_ol = ol_dir / 'sagittal' / f'{spec.specimen_id}.png'
        co_ol = ol_dir / 'coronal' / f'{spec.specimen_id}.png'

        reg_dir = spec.qc_grey_dirs
        ax_ol_r = natsorted((reg_dir / 'axial').iterdir(), key= lambda x: x.name)[-1]
        sa_ol_r = natsorted((reg_dir / 'sagittal').iterdir() ,key= lambda x: x.name)[-1]
        co_ol_r = natsorted((reg_dir / 'coronal').iterdir(), key= lambda x: x.name)[-1]

        f = plt.figure(figsize=(25, 12), dpi=96)

        show_(f, sa_ol, n_cols, 1)
        show_(f, co_ol, n_cols, 2)
        show_(f, ax_ol, n_cols, 3)
        show_(f, sa_ol_r, n_cols, 4, cmap='gray')
        show_(f, co_ol_r, n_cols, 5, cmap='gray')
        show_(f, ax_ol_r, n_cols, 6, cmap='gray', flip=True)
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