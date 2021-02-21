"""
Gicen a folder of stats results (each line in a subfolder), create
 heatmaps to summarize the hits across the specimens and at the line-level"""

from pathlib import Path
from typing import Iterable
import pandas as pd
from lama.stats.heatmap import heatmap
import matplotlib.pyplot as plt
from logzero import logger as logging


def heatmaps_for_permutation_stats(root_dir: Path):
    """
    This function works on the output of the premutation stats. For the non-permutation, may need to make a different
    function to deal with different directory layout
    """

    for line_dir in root_dir.iterdir():

        spec_dir = line_dir / 'specimen_level'
        spec_csvs = []

        for s_dir in spec_dir.iterdir():
            scsv = next(s_dir.iterdir())
            spec_csvs.append(scsv)

        try:
            line_hits_csv = next(line_dir.glob(f'{line_dir.name}_organ_volumes*csv'))
        except StopIteration:
            logging.error(f'cannot find stats results file in {str(line_dir)}')
            return

        line_specimen_hit_heatmap(line_hits_csv, spec_csvs, line_dir, line_dir.name)


def line_specimen_hit_heatmap(line_hits_csv: Path,
                              specimen_hits: Iterable[Path],
                              outdir: Path,
                              line: str,
                              sorter_csv=None):

    dfs = {}  # All line and speciemn hit dfs

    line_hits = pd.read_csv(line_hits_csv, index_col=0)
    dfs[line] = line_hits

    for spec_file in specimen_hits:
        d = pd.read_csv(spec_file, index_col=0)

        # d = d[d['significant_cal_p'] == True]

        # small_id = get_specimen_id(hits_file.name) for now use full name
        # dfs[small_id] = d
        dfs[spec_file.name] = d

    # get the superset of all hit labels
    hit_lables = set()
    for k, x in dfs.items():
        hit_lables.update(x[x['significant_cal_p'] == True].label_name)

    # For each hit table, keep only those in the hit superset and create heat_df
    t = []
    for line_or_spec, y in dfs.items():
        y = y[y['label_name'].isin(hit_lables)]
        y['label_num']= y.index
        y.set_index('label_name', inplace=True, drop=True)

        y.loc[y.significant_cal_p == False, 'mean_vol_ratio'] = None

        if 'mean_vol_ratio' in y:
            col_for_heatmap = 'mean_vol_ratio'
        else:
            col_for_heatmap = 'significant_cal_p'

        # Rename the column we are to display to the name of the specimen
        y.rename(columns={col_for_heatmap: line_or_spec}, inplace=True)
        t.append(y[[line_or_spec]])

    heat_df = pd.concat(t, axis=1)

    # if sorter_csv:
    #     # We have a csv with column abs(vol-diff) that we can use to sort results
    #     sorter = pd.read_csv(sorter_csv, index_col=0)
    #     s = 'abs(vol_diff)'
    #     heat_df = heat_df.merge(sorter[[s, 'label_name']], left_index=True, right_on='label_name', how='left')
    #     heat_df.sort_values(s, inplace=True, ascending=False)
    #     heat_df.set_index('label_name', drop=True, inplace=True)
    #     heat_df.drop(columns=[s], inplace=True)

    if len(heat_df) < 1:
        print(f"skipping {line}. No hits")
        return

    title = f'{line}\n  Line 5% FDR. Specimen 20% FDR'
    # Try to order lines by litter
    ids = list(heat_df.columns)
    line_id = ids.pop(0)

    ids.sort(key=lambda x: x[-2])
    sorted_ids = [line_id] + ids
    heat_df = heat_df[sorted_ids]

    # Shorten some of the longer organ names
    # label_meta = pd.read_csv('/mnt/bit_nfs/neil/Harwell_E14_5_latest/padded_target_ras/E14_5_atlas_v24_40_label_info_nouse.csv', index_col=0)
    # label_meta.set_index('label_name', inplace=True, drop=True)
    # heat_df = heat_df.merge(label_meta[['short_name']], how='left', left_index=True, right_index=True)

    def sn_mapper(x):
        sn = heat_df.loc[x, 'short_name']
        if not isinstance(sn, float): # float nan
            idx = sn
        else:
            idx = x
        return idx.lower()

    # heat_df.index = heat_df.index.map(sn_mapper)
    # heat_df.drop(columns='short_name', inplace=True)

    heatmap(heat_df, title=title, use_sns=True)

    plt.tight_layout()

    plt.savefig(outdir / f"{line}_organ_hit_heatmap.png")
    plt.close()


if __name__ == '__main__':
    spec_dir = Path('/mnt/bit_nfs/neil/impc_e15_5/phenotyping_tests/JAX_E15_5_test_120720/stats/archive/organ_vol_perm_091020/lines/Cox7c/specimen_level')
    spec_csvs = []

    for s_dir in spec_dir.iterdir():
        scsv = next(s_dir.iterdir())
        spec_csvs.append(scsv)
    line_specimen_hit_heatmap(Path('/mnt/bit_nfs/neil/impc_e15_5/phenotyping_tests/JAX_E15_5_test_120720/stats/archive/organ_vol_perm_091020/lines/Cox7c/Cox7c_organ_volumes_2020-10-09.csv'),
                              spec_csvs,
                              Path('/mnt/bit_nfs/neil/impc_e15_5/phenotyping_tests/JAX_E15_5_test_120720/stats/archive/organ_vol_perm_091020/lines/Cox7c'),
                        'Cox7c')