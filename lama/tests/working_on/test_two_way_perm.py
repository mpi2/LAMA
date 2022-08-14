# from lama.stats.permutation_stats.run_permutation_stats import run
import pytest
from pathlib import Path
from lama.scripts import lama_permutation_stats
from lama.stats.linear_model import lm_sm
from lama.stats.permutation_stats.run_permutation_stats import *
import lama.scripts.lama_permutation_stats

import random
from pathlib import Path
from datetime import date
import itertools
import pandas as pd
import numpy as np
from scipy.stats import zmap
from logzero import logger as logging
import yaml

from lama import common
from lama.stats.permutation_stats import distributions
from lama.stats.permutation_stats import p_thresholds
from lama.paths import specimen_iterator, get_specimen_dirs, LamaSpecimenData
from lama.qc.organ_vol_plots import make_plots, pvalue_dist_plots
from lama.common import write_array, read_array, init_logging, git_log, LamaDataException
from lama.stats.common import cohens_d
from lama.stats.penetrence_expressivity_plots import heatmaps_for_permutation_stats
from lama.stats.permutation_stats import run_permutation_stats

from lama.stats.permutation_stats.distributions import two_way_max_combinations, recursive_comb_maker, null_line
from lama.tests import (out_dir, wt_dir, mut_dir, treat_dir, inter_dir, label_meta, n_perm)

cfg = Path(
    "C:/Users/u5823099/Anaconda3/Lib/site-packages/lama/LAMA/lama/tests/configs/standard_stats/generate_data.toml")

stats_cfg = Path(
    "C:/Users/u5823099/Anaconda3/Lib/site-packages/lama/LAMA/lama/tests/configs/permutation_stats/perm_no_qc.yaml")


@pytest.mark.skip
def test_data_prep(two_way=True):
    """Just test if I can normalise stuff properly"""
    np.random.seed(999)
    init_logging(out_dir / 'stats.log')
    logging.info(git_log())
    logging.info(f'Running {__name__} with following commands\n{common.command_line_agrs()}')

    logging.info('Searching for staging data')
    wt_staging = get_staging_data(wt_dir)
    mut_staging = get_staging_data(mut_dir)

    logging.info('searching for organ volume data')
    wt_organ_vol = get_organ_volume_data(wt_dir)
    mut_organ_vol = get_organ_volume_data(mut_dir)
    if two_way:
        logging.info('Searching for two-way staging and organ volume data')
        treat_staging = get_staging_data(treat_dir)
        inter_staging = get_staging_data(inter_dir)
        treat_organ_vol = get_organ_volume_data(treat_dir)
        inter_organ_vol = get_organ_volume_data(inter_dir)
        two_way_data = [treat_staging, treat_organ_vol,
                        inter_staging, inter_organ_vol]

    data = prepare_data(wt_organ_vol,
                        wt_staging,
                        mut_organ_vol,
                        mut_staging,
                        label_meta=label_meta,
                        normalise_to_whole_embryo=True,
                        qc_file=None,
                        two_way=two_way,
                        two_way_data=two_way_data)

    # So pandas is a bit stupid as read/write of the file is causing type errors
    # How about we just look at
    print("Comparison")
    print(data)
    data.to_csv()
    #good_data = pd.read_csv('E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/input_data.csv', index_col=0)
    # check that the data is the same as a checked csv file
    #print(good_data)

    # print(data.astype(int).equals(data.astype(int)))


@pytest.mark.skip
def test_max_two_way_combinations():
    combs = two_way_max_combinations(num_wts=16)

    assert combs == 63063000


@pytest.mark.skip
def test_loop_iteration():
    lst = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    fin_results = []
    results = []

    for perm in range(10):
        full_combination = recursive_comb_maker(lst, 3, i=1, recurs_results=[])
        # print(full_combination)
        print(full_combination[1])
        fin_results.append(full_combination)
        # reset result list
        results = []
    print(fin_results[1])


@pytest.mark.skip
def test_generate_two_way_combinations():
    data = pd.read_csv('E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/input_data.csv', index_col=0)
    # can't really test this as of yet
    wt_indx_combinations = distributions.generate_random_two_way_combinations(data, 20)
    return wt_indx_combinations


@pytest.mark.skip
def test_lm_sm():
    data = pd.read_csv("C:/Users/u5823099/Anaconda3/Lib/site-packages/lama/LAMA/lama/tests/working_on/data.csv")
    info = pd.read_csv("C:/Users/u5823099/Anaconda3/Lib/site-packages/lama/LAMA/lama/tests/working_on/treat_info.csv")
    data = data.to_numpy()
    p, t = lm_sm(data=data, info=info, use_staging=True, two_way=True)
    print("p = ", p)
    print("t = ", t)


@pytest.mark.skip
def test_two_way_null_line():
    data = pd.read_csv('E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/input_data.csv', index_col=0)
    # can't really test this as of yet
    baselines = data[data['line'] == 'baseline']
    wt_indx_combinations = distributions.generate_random_two_way_combinations(data, 20)
    results = null_line(wt_indx_combinations, baselines, num_perms=20, two_way=True)

    print(type(results['x3']), results['x3'][0][0], type(results['x3'][0][0]))


@pytest.mark.skip
def test_two_way_null():
    data = pd.read_csv('E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/input_data.csv', index_col=0)

    line_null, specimen_null = distributions.null(input_data=data, num_perm=3, two_way=True)
    print(type(line_null['3'][0]), line_null['3'][0][0], type(line_null['3'][0][0]))
    print(type(specimen_null['22'][0]), specimen_null['22'][0][0], type(specimen_null['22'][0][0]))


@pytest.mark.skip
def test_two_way_alt():
    data = pd.read_csv('E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/input_data.csv', index_col=0)
    line_alt, spec_alt, line_alt_t, spec_alt_t = distributions.alternative(data, two_way=True)
    print(type(line_alt['3'][0]), line_alt['3'][0][0], type(line_alt['3'][0][0]))
    print(type(spec_alt['3'][0]), spec_alt['3'][0][0], type(spec_alt['3'][0][0]))


@pytest.mark.skip
def test_two_way_p_thresholds():
    """
     Testing the p_thresholds calculation
     -------
     TODO: Add more tests for different cases
     """

    def strip_x(dfs):
        for df in dfs:
            df.columns = [x.strip('x') for x in df.columns]

        return df

    data = pd.read_csv('E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/input_data.csv', index_col=0)

    baselines = data[data['line'] == 'baseline']
    wt_indx_combinations = distributions.generate_random_two_way_combinations(data, 100)
    results = null_line(wt_indx_combinations, baselines, num_perms=100, two_way=True)
    line_null = strip_x([results])

    # line_null, specimen_null = distributions.null(input_data=data, num_perm=10, two_way=True)

    line_alt, spec_alt, line_alt_t, spec_alt_t = distributions.alternative(data, two_way=True)

    thresh = p_thresholds.get_thresholds(line_null, line_alt, two_way=True)

    thresh.to_csv('E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/spec_out_threshs.csv')


@pytest.mark.skip
def test_two_spec_thresholds():
    two_way = True
    data = pd.read_csv('E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/input_data.csv', index_col=0)
    line_null, specimen_null = distributions.null(input_data=data, num_perm=3, two_way=True)

    line_alt, spec_alt, line_alt_t, spec_alt_t = distributions.alternative(data, two_way=True)

    # TODO: Don't hard-code this
    specimen_inter_nulls = specimen_null[specimen_null['3'].str.len() == 3]

    specimen_main_nulls = specimen_null[specimen_null['3'].str.len() == 1]
    specimen_geno_nulls, specimen_treat_nulls = np.vsplit(specimen_main_nulls, 2)

    specimen_inter_alt = spec_alt[spec_alt['3'].str.len() == 3]
    specimen_main_alt = spec_alt[spec_alt['3'].str.len() == 1]

    # TODO: Don't hard-code this

    specimen_geno_alt = specimen_main_alt[specimen_main_alt.index.str.contains("het")]
    specimen_treat_alt = specimen_main_alt[specimen_main_alt.index.str.contains("b6ku")]

    geno_thresholds = p_thresholds.get_thresholds(specimen_geno_nulls, specimen_geno_alt, two_way=two_way)
    treat_thresholds = p_thresholds.get_thresholds(specimen_treat_nulls, specimen_treat_alt, two_way=two_way)
    inter_thresholds = p_thresholds.get_thresholds(specimen_inter_nulls, specimen_inter_alt, two_way=two_way)

@pytest.mark.skip
def test_line_annotate():
    # Lines
    alt_file = Path('E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/distributions/specimen_geno_pvals.csv')
    thresholds_file = Path(
        'E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/distributions/specimen_geno_p_thresholds.csv')
    cond_dir = Path('E:/Bl6_data/211014_g_by_back')
    data = pd.read_csv('E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/input_data.csv', index_col=0)
    thresholds = pd.read_csv(thresholds_file, index_col=0)
    alt = pd.read_csv(alt_file, index_col=0)

    alt = alt.applymap(lambda x: np.array([float(i) for i in x.strip("[]").split()]) if "[" in x else x)


    #alt = alt.applymap(lambda x: np.array([float(i) for i in x.strip("[]").split()]))

    # run_permutation_stats.annotate(thresholds, alt, cond_dir, two_way=True, organ_volumes=data)

    # Specimens
    geno_hits = run_permutation_stats.annotate(thresholds, alt, cond_dir, is_line_level=False, two_way=False, organ_volumes=data, main_of_two_way=True)
    print(geno_hits)

@pytest.mark.skip
def test_add_significance():
    df = pd.read_csv('E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/test_df.csv', index_col=0)
    print(df)
    add_two_way_significance(df, 0.05)


@pytest.mark.skip
def test_two_way_plotting():
    data = pd.read_csv('E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/input_data.csv', index_col=0)

    label_info = Path('E:/Bl6_data/211014_g_by_back/target/E14_5_atlas_v24_43_label_info.csv')

    lines_root_dir = Path('E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output')

    normalise_to_whole_embryo = True
    voxel_size = 40

    data_for_plots = data.copy()
    data_for_plots.columns = [x.strip('x') for x in data_for_plots.columns]  # Strip any xs
    # If data has been normalised to WEV revert back for plots
    if normalise_to_whole_embryo:
        for col in data_for_plots.columns:
            if col.isdigit():
                data_for_plots[col] = data_for_plots[col] * data_for_plots['staging']

    make_plots(data_for_plots, label_info, lines_root_dir, voxel_size=voxel_size, two_way=True, skip_no_analysis=True)


@pytest.mark.skip
def test_dist_plots():
    line_null = pd.read_csv(
        "E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/distributions/null_line_dist_pvalues.csv")
    line_alt = pd.read_csv(
        "E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/distributions/alt_line_dist_pvalues.csv")
    line_organ_thresholds = pd.read_csv(
        "E:/Bl6_data/211014_g_by_back/permutation_stats/perm_output/distributions/line_organ_p_thresholds.csv")
    dist_plot_root = out_dir / 'distribution_plots'
    line_plot_dir = dist_plot_root / 'line_level'
    line_plot_dir.mkdir(parents=True, exist_ok=True)

    line_null_vals = pd.DataFrame([x[1:1000] for x in line_null.values])


    # line_null_vals.columns = line_null.columns
    line_null_vals.columns = line_null.drop(columns='Unnamed: 0').columns

    line_null_vals = line_null_vals.applymap(lambda x: np.array([float(i) for i in x.strip("[]").split()]))

    line_alt_vals = pd.DataFrame([x.strip("[]").split() for x in line_alt.values[0]])

    line_alt_vals = line_alt_vals.drop(line_alt_vals.index[0]).astype(float).transpose()

    line_alt_vals.columns = line_null_vals.columns

    pvalue_dist_plots(line_null_vals, line_alt_vals, line_organ_thresholds, line_plot_dir, label_meta_file=label_meta,
    two_way=True)


def test_two_way_heatmaps():
    label_info = Path('E:/Bl6_data/211014_g_by_back/target/E14_5_atlas_v24_43_label_info.csv')
    lines_root_dir = Path('E:/220607_two_way/permutation_stats/perm_output')
    heatmaps_for_permutation_stats(lines_root_dir, two_way=True,label_info_file=label_info)


@pytest.mark.skip
def test_permutation_stats():
    """
    Run the whole permutation based stats pipeline.
    Copy the output from a LAMA registrations test run, and increase or decrease the volume of the mutants so we get
    some hits

    """
    lama_permutation_stats.run(stats_cfg)