
from pathlib import Path
from lama.stats.standard_stats.radiomics import radiomics_job_runner
from lama.common import cfg_load
from lama.img_processing import normalise
import logging

def test_radiomics():

        c = cfg_load(Path("F:/220607_two_way/radiomics_output/generate_radiomics.toml"))

        target_dir = Path(c.get('target_dir'))

        labs_of_int = c.get('labs_of_int')

        norm_methods = [c.get('norm_methods')]


        norm_label = c.get('norm_label')

        spherify = c.get('spherify')

        ref_vol_path = Path(c.get('ref_vol_path'))

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
        radiomics_job_runner(target_dir, labs_of_int=labs_of_int,
                             normalisation_label=norm_label,
                             norm_method=norm_meths, spherify=spherify, ref_vol_path=ref_vol_path)