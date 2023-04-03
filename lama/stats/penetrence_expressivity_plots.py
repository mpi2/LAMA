"""
Gicen a folder of stats results (each line in a subfolder), create
 heatmaps to summarize the hits across the specimens and at the line-level"""

from pathlib import Path
from typing import Iterable
import pandas as pd
from lama.stats.heatmap import heatmap, clustermap
import matplotlib.pyplot as plt
from logzero import logger as logging
import numpy as np
from tqdm import tqdm

def heatmaps_for_permutation_stats(root_dir: Path, two_way: bool = False, label_info_file: Path = None, rad_plot: bool = False):
    """
    This function works on the output of the premutation stats. For the non-permutation, may need to make a different
    function to deal with different directory layout
    """
    # Yeah there should a be better way of me doing this but there is not
    if label_info_file:
        label_info = pd.read_csv(label_info_file, index_col=0)
        skip_analysis = True if 'no_analysis' in label_info.columns else False
        if skip_analysis:
            good_labels = label_info[label_info['no_analysis'] != True].label_name
    else:
        good_labels = None

    if two_way: # read data.csv  to get the conditions
        data_path = root_dir / "radiomics_data.csv" if rad_plot else root_dir / "input_data.csv"
        data = pd.read_csv(data_path, index_col=0)
        #if rad_plot:
        #    data.set_index(['vol'], inplace=True)
        group_info = data['line']

        # TODO: think whether to truly put mut_treat in main comparisons
        mut_names = data[group_info == 'mutants'].index
        treat_names = data[group_info == 'treatment'].index
        inter_names = data[group_info == 'mut_treat'].index

    for line_dir in root_dir.iterdir():
        # bodge way to fix it but she'll do
        if two_way:
            line_dir = root_dir / 'two_way'
            spec_dir = line_dir / 'specimen_level'

            file_lists = {
                'inter': [],
                'treat': [],
                'geno': []
            }

            for s_dir in spec_dir.iterdir():
                scsv = next(s_dir.iterdir())
                if s_dir.name in inter_names:
                    file_lists['inter'].append(scsv)
                elif s_dir.name in mut_names:
                    file_lists['geno'].append(scsv)
                elif s_dir.name in treat_names:
                    file_lists['treat'].append(scsv)

        else:
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

        if two_way:

            line_specimen_hit_heatmap(line_hits_csv, file_lists['geno'], line_dir, "geno", two_way=two_way,
                                      good_labels=good_labels, rad_plot=rad_plot)
            line_specimen_hit_heatmap(line_hits_csv, file_lists['treat'], line_dir, "treat", two_way=two_way,
                                      good_labels=good_labels, rad_plot=rad_plot)
            line_specimen_hit_heatmap(line_hits_csv, file_lists['inter'], line_dir, "inter", two_way=two_way,
                                      good_labels=good_labels, rad_plot=rad_plot)
            # there's only one iteration
            break
        else:
            line_specimen_hit_heatmap(line_hits_csv, spec_csvs, line_dir, line_dir.name, two_way=two_way,
                                    good_labels=good_labels)




def line_specimen_hit_heatmap(line_hits_csv: Path,
                              specimen_hits: Iterable[Path],
                              outdir: Path,
                              line: str,
                              sorter_csv=None,
                              two_way: bool = False, good_labels=None, rad_plot: bool = False):
    dfs = {}  # All line and speciemn hit dfs

    #line_hits = pd.read_csv(line_hits_csv, index_col=0)

    #dfs[line] = line_hits

    for spec_file in specimen_hits:
        d = pd.read_csv(spec_file, index_col=0)
        print("d", d)
        dfs[spec_file.name] = d

    # get the superset of all hit labels
    hit_lables = set()
    for k, x in tqdm(dfs.items()):

        # get significance c
        col = [_col for _col in x.columns if _col.__contains__("significant_cal")]
        # interaction has more than one p_val
        # [0] converts list to str
        col = "significant_cal_p_inter" if len(col) > 1 else col[0]

        if 'label_name' in x:
            # filtering for orgs of int
            if len(good_labels) > 1:

                good_hits = x[(x[col] == True) & (~x['no_analysis'].fillna(False))]

                good_hits = good_hits[good_hits['label_name'].isin(good_labels)].label_name

                hit_lables.update(good_hits)

            else:
                hit_lables.update(x[x[col] == True].label_name)
        else:

            hit_lables.update(x[x[col] == True].index.values)



    # For each hit table, keep only those in the hit superset and create heat_df
    t = []
    for line_or_spec, y in tqdm(dfs.items()):
        # If we have label_name, set as index. Otherwise leave label num as index
        if rad_plot:
            y = y[y['label_name'].isin(hit_lables)]
            # attach the label name to the index column
            #y['radiomic_name'] = [str(index.str.split("__")[0] + " " + row['label_name']).replace('_', ' ') for index, row in y.iterrows()]
            # thanks ChatGPT
            y['radiomic_name'] = y.apply(lambda row: f"{row.name.split('__')[0]} {row['label_name']}".replace('_', ' '),
                                         axis=1)
            #y.set_index('label', inplace=True, drop=True)
            y.set_index('radiomic_name', inplace=True, drop=True)

        elif 'label_name' in y:
            y = y[y['label_name'].isin(hit_lables)]
            y.set_index('label_name', inplace=True, drop=True)
            print("y", y)
        else:
            y.index = y.index.astype(str)

        y['label_num'] = y.index
        y.loc[y[col] == False, 'mean_vol_ratio'] = None

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
    ids.sort(key=lambda x: x[-3])
    sorted_ids = [line_id] + ids
    heat_df = heat_df[sorted_ids]

    if rad_plot:
        heat_df.columns = [col.rsplit('_', 3)[0] for col in heat_df.columns]

    try:
        if two_way:
            if not rad_plot:
                heat_df.columns = [x.split("org")[0] for x in heat_df.columns]

            if not heatmap(heat_df, title=title, use_sns=True, rad_plot=rad_plot):
                logging.info(f'Skipping heatmap for {line} as there are no results')

            plt.tight_layout()

            plt.savefig(outdir / f"{line}_organ_hit_heatmap.png")
            plt.close()

            #sns.clustermap needs non-nan values to calculate distances

            heat_df.dropna(how='all', inplace=True)
            heat_df.fillna(1, inplace=True)
            heat_df.clip(upper=2, inplace=True)

            if not clustermap(heat_df, title=title, use_sns=True, rad_plot=rad_plot):
                logging.info(f'Skipping heatmap for {line} as there are no results')

            plt.tight_layout()

            plt.savefig(outdir / f"{line}_organ_hit_clustermap.png")
            plt.close()


            logging.info("Creating Additional z-normalised plots")
            if not clustermap(heat_df, title=title, use_sns=True, rad_plot=rad_plot, z_norm=True):
                logging.info(f'Skipping heatmap for {line} as there are no results')

            plt.tight_layout()

            plt.savefig(outdir / f"{line}_organ_hit_clustermap_z_normed.png")
            plt.close()


        else:
            if not heatmap(heat_df, title=title, use_sns=True):
                logging.info(f'Skipping heatmap for {line} as there are no results')

            #plt.tight_layout()

            plt.savefig(outdir / f"{line}_organ_hit_heatmap.png")
            plt.close()
    except ValueError as e:
        print(e)
        logging.warn('No heatmap produced')


if __name__ == '__main__':
    spec_dir = Path(
        '/mnt/bit_nfs/neil/impc_e15_5/phenotyping_tests/JAX_E15_5_test_120720/stats/archive/organ_vol_perm_091020/lines/Cox7c/specimen_level')
    spec_csvs = []

    for s_dir in spec_dir.iterdir():
        scsv = next(s_dir.iterdir())
        spec_csvs.append(scsv)
    line_specimen_hit_heatmap(Path(
        '/mnt/bit_nfs/neil/impc_e15_5/phenotyping_tests/JAX_E15_5_test_120720/stats/archive/organ_vol_perm_091020/lines/Cox7c/Cox7c_organ_volumes_2020-10-09.csv'),
        spec_csvs,
        Path(
            '/mnt/bit_nfs/neil/impc_e15_5/phenotyping_tests/JAX_E15_5_test_120720/stats/archive/organ_vol_perm_091020/lines/Cox7c'),
        'Cox7c')
