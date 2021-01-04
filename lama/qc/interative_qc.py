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
from lama.qc.img_grid import HtmlGrid
import os

from lama.paths import get_specimen_dirs


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


def do_qc(root: Path, csv_or_dir: Path, html=False):
    """

    Parameters
    ----------
    root
        A Lama registrtion root for a line (or the baselines) containing multiple specimens
    csv
        Qc csv
    html
        if True, write to html for prettier formatiing
    """
    # Create series of images specimen-based view
    # mpl.rcParams['figure.raise_window']

    status_map = {
        'f': 'failed',
        'p': 'passed',
        'c': 'check'
    }

    if csv_or_dir.is_dir():
        interactive = False
        csv_or_dir = Path(csv_or_dir)
    else:
        interactive = True
        plt.ion()

        if csv_or_dir.is_file():
            df = pd.read_csv(csv_or_dir, index_col=0)
        else:
            df = pd.DataFrame(columns=['line', 'spec_id', 'status', 'comment'])

    # get first monitor. Could also look for largest monitor?
    mon = get_monitors()[0]
    width = mon.width
    width = width - (width * 0.1)  # Save som sapce for the console
    height = mon.height
    dpi = 100
    # What size does the figure need to be in inches to fit the image?
    figsize = width / float(dpi), height / float(dpi)
    f = plt.figure(figsize=figsize, dpi=dpi)

    grid = ImageGrid(f, 111,  # similar to subplot(111)
                     nrows_ncols=(3, 5),  # creates 2x2 grid of axes
                     axes_pad=0.1,  # pad between axes in inch.
                     share_all=False,
                     aspect=True
                     )

    for i, spec in enumerate(get_specimen_dirs(root)):

        try:
            spec.setup()
        except FileNotFoundError as e:
            print(f'Skipping {spec.specimen_id}\n{e}')
            continue

        if interactive:
            if spec.specimen_id in df.spec_id.values:
                print(f'skipping {spec.specimen_id}')
                continue
        print(f'QC {spec.specimen_id}')

        # Get the overlays. Use natsort to sort asceding by the slice number file suffix
        ol_dir = spec.qc_inverted_labels
        ax_ol = natsorted([x for x in (ol_dir / 'axial').iterdir()])
        sa_ol = natsorted([x for x in (ol_dir / 'sagittal').iterdir()])
        co_ol = natsorted([x for x in (ol_dir / 'coronal').iterdir()])

        reg_dir = spec.qc_grey_dirs
        ax_ol_r = natsorted((reg_dir / 'axial').iterdir(), key= lambda x: x.name)[-1]
        sa_ol_r = natsorted((reg_dir / 'sagittal').iterdir() ,key= lambda x: x.name)[-1]
        co_ol_r = natsorted((reg_dir / 'coronal').iterdir(), key= lambda x: x.name)[-1]

        ax_ol.append(ax_ol_r)
        sa_ol.append(sa_ol_r)
        co_ol.append(co_ol_r)

        all_p = chain(co_ol, ax_ol, sa_ol)

        i = 0

        if not html:  # mpl jpg output

            for ax, im in zip(grid, all_p):
                img = imread(im)
                if len(img.shape) == 2:
                    ax.imshow(img, cmap='gray')
                else:
                    ax.imshow(img)
                i += 1

            f.suptitle(f'{spec.line_id} - {spec.specimen_id}')
            f.tight_layout()

            if interactive:
                plt.show()
                # plt.pause(0.1)
                while True:
                    response = input("Enter status and optional comment: ")
                    status = response.split()[0]
                    if status not in status_map.keys():
                        print(f'status must be one of {status_map.keys()}')
                    else:
                        break

                comment = ' '.join(response.split()[1:])

                df.loc[len(df)] = [spec.line_id, spec.specimen_id, status, comment]
                df.to_csv(csv_or_dir)
            else:
                line_out_dir = csv_or_dir / spec.line_id
                line_out_dir.mkdir(exist_ok=True)
                spec_out_file = line_out_dir / f'{spec.specimen_id}.jpg'
                plt.savefig(spec_out_file, quality=80, optimize=True, progressive=True)

        else:  # Html output
            html_grid = HtmlGrid(spec.specimen_id)

            line_out_dir = csv_or_dir / spec.line_id
            line_out_dir.mkdir(exist_ok=True)
            spec_out_file = line_out_dir / f'{spec.specimen_id}.html'

            for g, img in enumerate(all_p):
                if g % 5 == 0:
                    html_grid.next_row('')
                img: Path
                img_rel = os.path.relpath(img, line_out_dir)
                html_grid.next_image(img_rel, '')

            html_grid.save(spec_out_file, width=400)


if __name__ =='__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser("Genate inverted label overaly plots. Multliple slices for each orietation")

    parser.add_argument('-i', '--indir', dest='indir',
                        help='A lama registration output directory containing one or more line directories',
                        required=True)
    parser.add_argument('-c', '--csv', dest='csv',
                        help='Use script interatively. View one specimen at a time and write qc status to this csv',
                        required=False, default=None)
    parser.add_argument('-d', '--dir', dest='dir',
                        help='Use script non-interatively. Write qc images to dir in line subfolders',
                        required=False, default=None)
    parser.add_argument('--html', dest='html',
                        help='Use with -d. write to html instead of using matplotlib for larger images',
                        required=False, default=None, action='store_true')
    parser.add_argument('-s', '--single_file', dest='single_file',
                        help='Embed images within the the html file (Use --html option to). Makes portable html files'
                             'that are not tied to the LAMA folder',
                        required=False, default=None, action='store_true')
    args = parser.parse_args()

    if args.csv is None and args.dir is None:
        sys.exit("Supply either a csv path (-c) or a directory (-d)")
    if args.csv is not None and args.dir is not None:
        sys.exit("Supply either a csv path (-c) or a directory (-d)")

    if args.csv is not None:
        csv_or_dir = Path(args.csv)
    else:
        csv_or_dir = Path(args.dir)
        if not csv_or_dir.is_dir():
            sys.exit('-d must be a directory')

    do_qc(Path(args.indir), csv_or_dir, args.html)