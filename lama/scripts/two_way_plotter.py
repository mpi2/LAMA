from lama.utilities import combine_spec_csv
from lama import common
import logging
import subprocess as sub
from pathlib import Path
PLOT_SCRIPT = str(common.lama_root_dir / 'stats' / 'rscripts' / 'two_way_plot.R')

def main(root_dir):
    # combine specimen organ volumes / staging volumes
    logging.info("Plotting Two-way Standard Stats")
    combine_spec_csv.main(root_dir)

    # run plotting script
    cmd = ['Rscript',
           PLOT_SCRIPT,
           str(root_dir / "full_organs.csv"),
           str(root_dir / "full_staging.csv"),
           str(root_dir.parent / 'target' / 'E14_5_atlas_v24_43_label_info'),
           40.0
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