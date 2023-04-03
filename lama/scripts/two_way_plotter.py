from lama.utilities import combine_spec_csv, extract_label, extract_registrations
from lama import common
import logging
import subprocess as sub
from pathlib import Path
import os
import tempfile
PLOT_SCRIPT = str(common.lama_root_dir / 'stats' / 'rscripts' / 'two_way_plot.R')


def main():
    import argparse

    parser = argparse.ArgumentParser("plot way data and get segmentations of interest")
    parser.add_argument('-i', dest='root_dir',
                        help='Folder with Registration results, (i.e. wild_type_and_mutant_data)',
                        required=True)
    parser.add_argument('-l', dest='labs',
                        help='labels of interest (list of numbers - look at label_info for numbers)',
                        required=True)

    args = parser.parse_args()
    root_dir = Path(args.root_dir)

    labs = [float(i) for i in args.labs.split(",")] if "," in args.labs else float(args.labs)

    logging.info("Extracting Volumes and Labels of Interest")
    extract_registrations.main(root_dir)
    extract_label.main(root_dir, labs)


    logging.info("Plotting Two-way Standard Stats")
    # combine specimen organ volumes / staging volumes
    combine_spec_csv.main(root_dir)

    organ_vol_file = tempfile.NamedTemporaryFile().name
    staging_file = tempfile.NamedTemporaryFile().name
    label_info_file = tempfile.NamedTemporaryFile().name
    groups_file = tempfile.NamedTemporaryFile().name



    # run plotting script

    print(os.path.exists(PLOT_SCRIPT))
    print(str(root_dir / "full_organs.csv"))
    print(str(root_dir / "full_staging.csv"))
    print(str(root_dir.parent / 'target' / 'E14_5_atlas_v24_43_label_info.csv'))
    cmd = ['Rscript',
           PLOT_SCRIPT,
           str(root_dir / "full_organs.csv"),
           str(root_dir / "full_staging.csv"),
           str(root_dir.parent / 'target' / 'E14_5_atlas_v24_43_label_info.csv'),
           '40.0'
           ]
    try:
        sub.check_output(cmd)
        logging.info('R plotting suceeded')
    except sub.CalledProcessError as e:
        msg = "R plotting failed {}".format(e)
        logging.exception(msg)
        raise RuntimeError(msg)


if __name__ == '__main__':
    main()
