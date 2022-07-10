
from lama.stats.standard_stats.radiomics import radiomics_job_runner
import pandas as pd
import logging
import pathlib as Path
from lama.common import cfg_load
from lama.img_processing import normalise

def main():
    import argparse

    parser = argparse.ArgumentParser("Schedule LAMA jobs")
    #TODO: make a config
    parser.add_argument('-c', '--config', dest='config', help='lama.yaml config file',
                        required=True)
    parser.add_argument('-r', '--root_dir', dest='root_dir', help='The root directory containing the input folders',
                        required=True)

    parser.add_argument('-m', '--make_job_file', dest='make_job_file', help='Run with this option forst to crate a job file',
                    action='store_true', default=False)


    args = parser.parse_args()

    try:
        # lets just get the config here - it's not that big right now

        c = cfg_load(Path(args.config))

        target_dir = c.get('target_dir')

        labs_of_int = c.get('labs_of_int')

        norm_methods = c.get('norm_methods')

        norm_label = c.get('norm_label')

        spherify = c.get('spherify')

        norm_dict = {
            "histogram": normalise.IntensityHistogramMatch(),
            "N4": normalise.IntensityN4Normalise(),
            "subtraction": normalise.NonRegMaskNormalise()
        }

        try:
            norm_meths = [norm_dict[x] for x in norm_methods]
        except KeyError:
            norm_meths = None

        radiomics_job_runner(target_dir, labs_of_int=labs_of_int,
                             normalisation_label=norm_label,
                             norm_method=norm_meths, spherify=spherify)
    except pd.errors.EmptyDataError as e:
        logging.exception(f'pandas read failure {e}')


if __name__ == '__main__':
    main()