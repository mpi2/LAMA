from lama.lama_radiomics.radiomics import radiomics_job_runner
import pandas as pd
import logging
from pathlib import Path
from lama.common import cfg_load
from lama.img_processing import normalise


def main():
    import argparse

    parser = argparse.ArgumentParser("Schedule LAMA jobs")

    parser.add_argument('-c', '--config', dest='config', help='lama.yaml config file',
                        required=True)

    parser.add_argument('-m', '--make_job_file', dest='make_job_file',
                        help='Run with this option forst to crate a job file',
                        action='store_true', default=False)

    args = parser.parse_args()

    try:
        # lets just get the config here - it's not that big right now

        c = cfg_load(Path(args.config))
        # c = cfg_load(Path("E:/220607_two_way/radiomics_output/generate_radiomics.toml"))

        target_dir = Path(c.get('target_dir'))

        labs_of_int = c.get('labs_of_int')

        norm_methods = [c.get('norm_methods')]

        print("from c.get", norm_methods)

        norm_label = c.get('norm_label')

        spherify = c.get('spherify')

        fold = c.get('fold')

        ref_vol_path = Path(c.get('ref_vol_path')) if c.get('ref_vol_path') is not None else None

        norm_dict = {
            "histogram": normalise.IntensityHistogramMatch(),
            "N4": normalise.IntensityN4Normalise(),
            "subtraction": normalise.NonRegMaskNormalise(),
            "none": None
        }
        try:
            norm_meths = [norm_dict[str(x)] for x in norm_methods]



        except KeyError as e:
            print(e)

            norm_meths = None
        logging.info("Starting Radiomics")

        print(norm_meths)
        radiomics_job_runner(target_dir, labs_of_int=labs_of_int, norm_method=norm_meths, spherify=spherify,
                             ref_vol_path=ref_vol_path, norm_label=norm_label, make_job_file=args.make_job_file, fold=fold)
    except pd.errors.EmptyDataError as e:
        logging.exception(f'pandas read failure {e}')


if __name__ == '__main__':
    main()
