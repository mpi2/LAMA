"""
09/04/19
data_loader in vpv.util can now load in a a series of volumes and setup the views in a predifned way.

This script will generate a config for each line.

Example output file
-------------------

labels_dir = 'inverted_labels/similarity'
vol_dir = 'registrations/rigid'
orientation = 'sagittal'


[top]
specimens = ['/mnt/IMPC_research/neil/E14.5/baselines/output/baseline/20150916_RBMS1_E14.5_14.3f_WT_XY_rec_scaled_4.6878_pixel_14',
'/mnt/IMPC_research/neil/E14.5/baselines/output/baseline/20170214_1200014J11RIK_E14.5_1.5h_WT_XX_REC_scaled_4.7297_pixel_13',
'/mnt/IMPC_research/neil/E14.5/baselines/output/baseline/20140121_RIC8B_E14.5_15.4b_wt_xy_rec_scaled_3.125_pixel_14']

[bottom]
specimens = ['/mnt/IMPC_research/neil/E14.5/mutants/output/1200014J11RIK/20170214_1200014J11RIK_E14.5_1.5f_HOM_XX_REC_scaled_4.7297_pixel_13.9999',
'/mnt/IMPC_research/neil/E14.5/mutants/output/1200014J11RIK/20170214_1200014J11RIK_E14.5_2.4c_HOM_XX_REC_scaled_4.7297_pixel_13.9999',
'/mnt/IMPC_research/neil/E14.5/mutants/output/1200014J11RIK/20170214_1200014J11RIK_E14.5_2.4i_HOM_XY_REC_scaled_4.7297_pixel_13.9999']


Using output file to load volumes in VPV
----------------------------------------
todo


"""

from pathlib import Path
import toml
import pandas as pd


def write_config(wt_paths, mut_paths, out_path):
    """

    Parameters
    ----------
    """
    config = dict(
    labels_dir = 'inverted_labels/similarity',
    vol_dir = 'registrations/rigid',
    orientation = 'sagittal')

    config['top'] = {'specimens': [str(x) for x in wt_paths]}
    config['bottom'] = {'specimens': [str(x) for x in mut_paths]}

    with open(out_path, 'w') as fh:
        fh.write(toml.dumps(config))


def run(wt_root, mut_root, wt_staging, out_dir):

    wt_staging = pd.read_csv(wt_staging, index_col=0)

    for line_stats_dir in mut_root.iterdir():

        print(line_stats_dir.name)

        wt_specimen_dirs = []
        mut_specimen_dirs = []

        wt_staging['used'] = False

        line_reg_dir = mut_root / line_stats_dir.name

        # Get each specimen rigid image. Get nearest staged available wt image
        if not line_reg_dir.is_dir():
            continue

        for specimen_dir in line_reg_dir.iterdir():
            if str(specimen_dir).endswith('_'):
                continue # Legacy stats were put there

            mut_specimen_dirs.append(specimen_dir)
            spec_stage_file = specimen_dir /'output' / 'staging_info_volume.csv'
            spec_stage_df = pd.read_csv(spec_stage_file, index_col=0)
            spec_stage = spec_stage_df.value.values[0]

            # Get a stage-matched wt
            wt_to_use = abs(wt_staging[wt_staging.used != True][['staging']] - spec_stage).sort_values(by='staging').iloc[0].name
            wt_staging.at[wt_to_use, 'used'] = True

            wt_dir = wt_root / wt_to_use
            wt_specimen_dirs.append(wt_dir)

        out_path = out_dir / f'{line_stats_dir.name}.toml'
        write_config(wt_specimen_dirs, mut_specimen_dirs, out_path)


if __name__ == '__main__':
    wt_root = Path('/mnt/IMPC_research/neil/E14.5/baselines/output/baseline')
    mut_root = Path('/mnt/IMPC_research/neil/E14.5/mutants/output')
    wt_staging = Path('/mnt/IMPC_research/neil/E14.5/baselines/output/staging_info_volume.csv')
    out_dir = Path('/mnt/IMPC_research/neil/E14.5/img_loaders')



    run(wt_root, mut_root, wt_staging, out_dir)










