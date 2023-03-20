
from pathlib import Path

import joblib

from lama.lama_radiomics.radiomics import radiomics_job_runner
from lama import common
import os
from lama.common import cfg_load
from lama.img_processing import normalise
from logzero import logger as logging
import seaborn as sns
import pandas as pd
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pytest
import numpy as np
from lama.lama_radiomics import feature_reduction
from lama.scripts import lama_machine_learning
import pacmap
from lama.scripts import lama_permutation_stats
from lama.lama_radiomics import radiomics

import SimpleITK as sitk
stats_cfg = Path(
    "C:/LAMA/lama/tests/configs/permutation_stats/perm_no_qc.yaml")

from lama.stats.permutation_stats.run_permutation_stats import get_radiomics_data

def test_denoising():
    file_path = Path("E:/220204_BQ_dataset/221218_BQ_run/registrations/rigid/flipped/200721_MPTLVo3_CT_4T1_Ms_D7_C1_002.nrrd")

    img = common.LoadImage(file_path).img
    f_out = "E:/220204_BQ_dataset/221218_BQ_run/test.nrrd"

    result = radiomics.denoise(img)
    sitk.WriteImage(result, f_out)





def test_radiomics():
        cpath =  Path('C:/LAMA/lama/tests/configs/lama_radiomics/radiomics_config_gina.toml')
        c = cfg_load(cpath)

        target_dir = Path(c.get('target_dir'))

        labs_of_int = c.get('labs_of_int')

        norm_methods = c.get('norm_methods')


        norm_label = c.get('norm_label')

        spherify = c.get('spherify')

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
        radiomics_job_runner(target_dir, labs_of_int=labs_of_int, norm_method=normalise.IntensityHistogramMatch(), norm_label=norm_label,
                             spherify=spherify, ref_vol_path=ref_vol_path, make_job_file=False)


def test_radiomic_plotting():
    _dir = Path("E:/220607_two_way/g_by_back_data/radiomics_output/features/")

    file_names = [spec for spec in common.get_file_paths(folder=_dir, extension_tuple=".csv")]
    file_names.sort()

    data = [pd.read_csv(spec, index_col=0).dropna(axis=1) for spec in file_names]

    abnormal_embs = ['22300_e8','22300_e6', '50_e5']

    for i, df in enumerate(data):
        df.index.name = 'org'
        df.name = str(file_names[i]).split(".")[0].split("/")[-1]
        df['genotype'] = 'HET' if 'het' in str(file_names[i]) else 'WT'
        df['background'] = 'C57BL6N' if (('b6ku' in str(file_names[i]))|('BL6' in str(file_names[i]))) else \
            'F1' if ('F1' in str(file_names[i])) else 'C3HHEH'

        df['HPE'] = 'abnormal' if any(map(str(file_names[i]).__contains__, abnormal_embs)) else 'normal'

    data = pd.concat(
        data,
        ignore_index=False, keys=[os.path.splitext(os.path.basename(spec))[0] for spec in file_names],
        names=['specimen', 'org'])

    line_file = _dir.parent / "full_results.csv"

    org_dir =_dir.parent / "organs"

    os.makedirs(org_dir, exist_ok=True)
    print(data.columns)

    for org in data.index.get_level_values('org').unique():
        data[data.index.get_level_values('org') == org].to_csv(str(org_dir)+"/results_" + str(org)+ ".csv")

    data.to_csv(line_file)

    data_subset = data.select_dtypes(include=np.number)

    data_subset = data_subset.apply(lambda x: (x - x.mean()) / x.std(), axis=0)
    data_subset = data_subset.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

    embedding = pacmap.PaCMAP(n_components=2, n_neighbors=10, MN_ratio=0.5, FP_ratio=2.0, num_iters=20000, verbose=1)



    #print(data_subset.dropna(axis='columns'))

    results = embedding.fit_transform(data_subset.dropna(axis='columns'))

    color_class = data.index.get_level_values('org')

    # fig, ax = plt.subplots(figsize=[55, 60])
    # cluster.tsneplot(score=tsne_results, show=True, theme='dark', colorlist=color_class)

    data['PaCMAP-2d-one'] = results[:, 0]
    data['PaCMAP-2d-two'] = results[:, 1]
    data['org'] = data.index.get_level_values('org')
    data['specimen'] = data.index.get_level_values('specimen')
    data['condition'] = data['genotype'] + "_" + data['background']


    fig, ax = plt.subplots(figsize=[56, 60])
    #data = data[data['condition'] == 'WT_C3HHEH']

    g = sns.lmplot(
        x="PaCMAP-2d-one", y="PaCMAP-2d-two",
        data=data,
        #col_order=['normal', 'abnormal'],
        col='condition',
        col_wrap=2,
        hue="org",
        palette='husl',
        fit_reg=False)
    g.set(ylim=(np.min(data['PaCMAP-2d-two'])-10, np.max(data['PaCMAP-2d-two'])+10),
          xlim=(np.min(data['PaCMAP-2d-one'])-10, np.max(data['PaCMAP-2d-one'])+10))


    plt.savefig("E:/220607_two_way/g_by_back_data/radiomics_output/features/radiomics_2D_PaCMAP_all_cond_v2.png")
    plt.close()

    fig, ax = plt.subplots(figsize=[56, 60])
    wt_c3h_data = data[data['condition'] == 'WT_C3HHEH']

    g = sns.lmplot(
        x="PaCMAP-2d-one", y="PaCMAP-2d-two",
        data=wt_c3h_data,
        # col_order=['normal', 'abnormal'],
        col='org',
        col_wrap=5,
        hue="org",
        palette='husl',
        fit_reg=False)
    g.set(ylim=(np.min(data['PaCMAP-2d-two']) - 10, np.max(data['PaCMAP-2d-two']) + 10),
          xlim=(np.min(data['PaCMAP-2d-one']) - 10, np.max(data['PaCMAP-2d-one']) + 10))

    plt.savefig("E:/220607_two_way/g_by_back_data/radiomics_output/features/radiomics_2D_PaCMAP_C3H_wt_org_v2.png")
    plt.close()

    fig, ax = plt.subplots(figsize=[56, 60])
    g = sns.lmplot(
        x="PaCMAP-2d-one", y="PaCMAP-2d-two",
        data=wt_c3h_data,
        # col_order=['normal', 'abnormal'],
        #col='specimen',
        #col_wrap=5,
        hue="org",
        palette='husl',
        fit_reg=False)
    g.set(ylim=(np.min(data['PaCMAP-2d-two']) - 10, np.max(data['PaCMAP-2d-two']) + 10),
          xlim=(np.min(data['PaCMAP-2d-one']) - 10, np.max(data['PaCMAP-2d-one']) + 10))

    plt.savefig("E:/220607_two_way/g_by_back_data/radiomics_output/features/radiomics_2D_PaCMAP_C3H_map_v2.png")
    plt.close()

    het_c3h_data = data[data['condition'] == 'HET_C3HHEH']
    fig, ax = plt.subplots(figsize=[56, 60])
    g = sns.lmplot(
        x="PaCMAP-2d-one", y="PaCMAP-2d-two",
        data=het_c3h_data,
        # col_order=['normal', 'abnormal'],
        col='specimen',
        col_wrap=5,
        hue="org",
        palette='husl',
        fit_reg=False)
    g.set(ylim=(np.min(data['PaCMAP-2d-two']) - 10, np.max(data['PaCMAP-2d-two']) + 10),
          xlim=(np.min(data['PaCMAP-2d-one']) - 10, np.max(data['PaCMAP-2d-one']) + 10))

    plt.savefig("E:/220607_two_way/g_by_back_data/radiomics_output/features/radiomics_2D_PaCMAP_C3H_hets_v2.png")
    plt.close()

    wt_b6_data = data[data['condition'] == 'WT_C57BL6N']
    fig, ax = plt.subplots(figsize=[56, 60])
    g = sns.lmplot(
        x="PaCMAP-2d-one", y="PaCMAP-2d-two",
        data=wt_b6_data,
        # col_order=['normal', 'abnormal'],
        #col='specimen',
        #col_wrap=5,
        hue="org",
        palette='husl',
        fit_reg=False)
    g.set(ylim=(np.min(data['PaCMAP-2d-two']) - 10, np.max(data['PaCMAP-2d-two']) + 10),
          xlim=(np.min(data['PaCMAP-2d-one']) - 10, np.max(data['PaCMAP-2d-one']) + 10))

    plt.savefig("E:/220607_two_way/g_by_back_data/radiomics_output/features/radiomics_2D_PaCMAP_b6_map_v2.png")
    plt.close()

    wt_f1_data = data[data['condition'] == 'WT_C57BL6N']
    fig, ax = plt.subplots(figsize=[56, 60])
    g = sns.lmplot(
        x="PaCMAP-2d-one", y="PaCMAP-2d-two",
        data=wt_f1_data,
        # col_order=['normal', 'abnormal'],
        # col='specimen',
        # col_wrap=5,
        hue="org",
        palette='husl',
        fit_reg=False)
    g.set(ylim=(np.min(data['PaCMAP-2d-two']) - 10, np.max(data['PaCMAP-2d-two']) + 10),
          xlim=(np.min(data['PaCMAP-2d-one']) - 10, np.max(data['PaCMAP-2d-one']) + 10))

    plt.savefig("E:/220607_two_way/g_by_back_data/radiomics_output/features/radiomics_2D_PaCMAP_f1_map_v2.png")
    plt.close()

    fig, ax = plt.subplots(figsize=[56, 60])
    g = sns.lmplot(
        x="PaCMAP-2d-one", y="PaCMAP-2d-two",
        data=wt_b6_data,
        # col_order=['normal', 'abnormal'],
        col='specimen',
        col_wrap=5,
        hue="org",
        palette='husl',
        fit_reg=False)
    g.set(ylim=(np.min(data['PaCMAP-2d-two']) - 10, np.max(data['PaCMAP-2d-two']) + 10),
          xlim=(np.min(data['PaCMAP-2d-one']) - 10, np.max(data['PaCMAP-2d-one']) + 10))

    plt.savefig("E:/220607_two_way/g_by_back_data/radiomics_output/features/radiomics_2D_PaCMAP_b6_specs_v2.png")
    plt.close()

    het_b6_data = data[data['condition'] == 'HET_C57BL6N']
    fig, ax = plt.subplots(figsize=[56, 60])
    g = sns.lmplot(
        x="PaCMAP-2d-one", y="PaCMAP-2d-two",
        data=het_b6_data,
        # col_order=['normal', 'abnormal'],
        col='specimen',
        col_wrap=5,
        hue="org",
        palette='husl',
        fit_reg=False)
    g.set(ylim=(np.min(data['PaCMAP-2d-two']) - 10, np.max(data['PaCMAP-2d-two']) + 10),
          xlim=(np.min(data['PaCMAP-2d-one']) - 10, np.max(data['PaCMAP-2d-one']) + 10))

    plt.savefig("E:/220607_two_way/g_by_back_data/radiomics_output/features/radiomics_2D_PaCMAP_b6_hets_v2.png")
    plt.close()

    het_f1_data = data[data['condition'] == 'F1']
    fig, ax = plt.subplots(figsize=[56, 60])
    g = sns.lmplot(
        x="PaCMAP-2d-one", y="PaCMAP-2d-two",
        data=het_f1_data,
        # col_order=['normal', 'abnormal'],
        col='specimen',
        col_wrap=5,
        hue="org",
        palette='husl',
        fit_reg=False)
    g.set(ylim=(np.min(data['PaCMAP-2d-two']) - 10, np.max(data['PaCMAP-2d-two']) + 10),
          xlim=(np.min(data['PaCMAP-2d-one']) - 10, np.max(data['PaCMAP-2d-one']) + 10))

    plt.savefig("E:/220607_two_way/g_by_back_data/radiomics_output/features/radiomics_2D_PaCMAP_f1_hets_v2.png")
    plt.close()

    fig, ax = plt.subplots(figsize=[56, 60])
    g = sns.lmplot(
        x="PaCMAP-2d-one", y="PaCMAP-2d-two",
        data=data,
        col_order=['normal', 'abnormal'],
        col='HPE',
        row='condition',
        hue="org",
        palette='husl',
        fit_reg=False)
    g.set(ylim=(np.min(data['PaCMAP-2d-two']) - 10, np.max(data['PaCMAP-2d-two']) + 10),
          xlim=(np.min(data['PaCMAP-2d-one']) - 10, np.max(data['PaCMAP-2d-one']) + 10))

    plt.savefig("E:/220607_two_way/g_by_back_data/radiomics_output/features/radiomics_2D_PaCMAP_HPE_v2.png")
    plt.close()



    #fig, ax = plt.subplots((len(data.index.levels[-1]) // 5) + 1, 5, figsize=[60, 80], sharex=True, sharey=True)

    #for i, row in enumerate(data.index.levels[-1]):
    #    fig, ax = plt.subplots(figsize=[27, 30])
    #    print(row)
    #    sns.lmplot(
    #        x="tsne-3d-one", y="tsne-2d-two",
    #        hue="org",
    #        col="condition",
    #        palette="husl",
    #        data=data.loc[row],
    #        legend="full",
    #        fit_reg=False,
    #        legend_out=True)
    #    file_name = "E:/220606_two_way/g_by_back_data/radiomics_output/sample_features/" + str(row) + ".png"
    #    plt.savefig(file_name)
    #    plt.close()

    #    sns.relplot(
    #        x="tsne-3d-one", y="tsne-2d-two",
    #        data=data,
    #        hue="org",
    #        palette='husl',
    #        alpha=-1.3,
    #        ax=ax[(i+0)//5, (i+1)%5])


    #fig.savefig("E:/220720_Amrit_radiomics/radiomics_2D_tsne_overlay.png")
    #plt.close()

    # remove diagnostics
    # data.index = data['specimen']
    # print(data.index.str.rsplit('_', 1))
    # data = data[data.columns.drop(list(data.filter(regex="diagnostics")))]

    # _metadata = pd.DataFrame(data.index.str.rsplit('_', 1))

    # print(_metadata)

    # _metadata[['Embryo','Genotype']] = pd.DataFrame(_metadata.specimen.tolist(), index=_metadata.index)

    # print(_metadata)

    # _metadata = _metadata.drop(columns=['specimen'])

    # _metadata.reset_index(inplace=True, drop=True)
    # data.reset_index(inplace=True, drop=True)

    # data=data.drop(columns=['specimen'])

    # print(data)
    # umap_organs(data, Path("E:/Bl5_data/211014_g_by_back/umap.png"), _metadata=_metadata)


    data.drop(['org', 'condition', 'HPE'], axis=1, inplace=True)

    data = data.select_dtypes(include=np.number)

    #data = data.apply(lambda x: (x - x.mean()) / x.std(), axis=0)

    data = data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)


    #data = data.apply(lambda x: (x - x.mean()) / x.std(), axis=0)

    #sns.set(font_scale=0.5)

    print("Data after drop", data)

    print(data.index.get_level_values('org'), len(data.index.get_level_values('org')))

    #data = data[~np.isin(data.index.get_level_values('org'), 27.0)]

    print(data)
    data.drop(27.0, level=1, axis=0, inplace=True)
    print(data, any(np.isin(data.index.get_level_values('org'), 27.0)))
    print(data.index.get_level_values('org'))

    for i, org in enumerate(data.index.levels[1]):
        fig, ax = plt.subplots(figsize=[14, 15])
        #sns.set(font_scale=0.5)
        o_data = data[np.isin(data.index.get_level_values('org'), org)]
        o_data.dropna(axis=1, inplace=True)

        if org == 27.0:
            continue

        print(org)

        o_data.drop(o_data.std()[(o_data.std() == 0)].index, axis=1, inplace=True)


        #o_data.to_csv("E:/org_df.csv")

        #import scipy.spatial.distance as ssd


        sns.clustermap(o_data.T,
                    figsize=[148.5, 210],
                    dendrogram_ratio=0.1,
                    colors_ratio=0.1,
                    #z_score=0,
                    metric="correlation",
                    #cmap=sns.diverging_palette(250, 15, l=70, s=400, sep=40, n=512, center="light", as_cmap=True),
                    #cbar_kws={'Genotype': 'Background'},
                    square=True,
                    xticklabels=True,
                    yticklabels=True)
        plt.tight_layout()
        plt.savefig("E:/220607_two_way/g_by_back_data/radiomics_output/wt_C3H_heatmap_"+str(org)+".png")
        plt.close()


def test_BQ_concat():
    _dir = Path("Z:/jcsmr/ROLab/Experimental data/Radiomics/Workflow design and trial results/Kyle Drover analysis/220617_BQ_norm_stage_full/sub/sub_normed_features.csv")
    #_dir = Path("E:/220913_BQ_tsphere/inputs/features/")

    # file_names = [spec for spec in common.get_file_paths(folder=_dir, extension_tuple=".csv")]
    # file_names.sort()
    #
    # data = [pd.read_csv(spec, index_col=0).dropna(axis=1) for spec in file_names]
    #
    # data = pd.concat(
    #     data,
    #     ignore_index=False, keys=[os.path.splitext(os.path.basename(spec))[0] for spec in file_names],
    #     names=['specimen', 'label'])
    #
    # data['specimen'] = data.index.get_level_values('specimen')
    #
    # _metadata = data['specimen'].str.split('_', expand=True)
    #
    #
    #
    # _metadata.columns = ['Date', 'Exp', 'Contour_Method', 'Tumour_Model', 'Position', 'Age',
    #                                                                         'Cage_No.', 'Animal_No.']
    #
    #
    #
    #
    # _metadata.reset_index(inplace=True, drop=True)
    # data.reset_index(inplace=True, drop=True)
    # features = pd.concat([_metadata, data], axis=1)
    #
    # features.index.name = 'scanID'
    #
    # print(features)
    #
    # print(str(_dir.parent / "full_results.csv"))
    #
    # features.to_csv(str(_dir.parent / "full_results.csv"))

    features = pd.read_csv(_dir)
    features = features[features.columns.drop(list(features.filter(regex="diagnostics")))]
    features.drop(["scanID"], axis=1, inplace=True)
    feature_reduction.main(features, org = None, rad_file_path = Path(_dir.parent / "full_results.csv"))

def test_BQ_mach_learn():
    _dir = Path("C:/test/features/")

    file_names = [spec for spec in common.get_file_paths(folder=_dir, extension_tuple=".csv")]
    file_names.sort()

    data = [pd.read_csv(spec, index_col=0).dropna(axis=1) for spec in file_names]

    data = pd.concat(
        data,
        ignore_index=False, keys=[os.path.splitext(os.path.basename(spec))[0] for spec in file_names],
        names=['specimen', 'label'])

    data['specimen'] = data.index.get_level_values('specimen')

    _metadata = data['specimen'].str.split('_', expand=True)



    _metadata.columns = ['Date', 'Exp', 'Contour_Method', 'Tumour_Model', 'Position', 'Age',
                                                                            'Cage_No.', 'Animal_No.']




    _metadata.reset_index(inplace=True, drop=True)
    data.reset_index(inplace=True, drop=True)
    features = pd.concat([_metadata, data], axis=1)

    features.index.name = 'scanID'

    print(features)

    print(str(_dir.parent / "full_results.csv"))

    features.to_csv(str(_dir.parent / "full_results.csv"))

    feature_reduction.main(features, org = None, rad_file_path = Path(_dir.parent / "full_results.csv"))


def test_BQ_mach_learn_non_tum():
    _dir = Path("E:/220919_non_tum/features/")

    file_names = [spec for spec in common.get_file_paths(folder=_dir, extension_tuple=".csv")]
    file_names.sort()

    data = [pd.read_csv(spec, index_col=0).dropna(axis=1) for spec in file_names]

    data = pd.concat(
        data,
        ignore_index=False, keys=[os.path.splitext(os.path.basename(spec))[0] for spec in file_names],
        names=['specimen', 'label'])

    data['specimen'] = data.index.get_level_values('specimen')

    _metadata = data['specimen'].str.split('_', expand=True)



    _metadata.columns = ['Date', 'Exp', 'Contour_Method', 'Tumour_Model', 'Position', 'Age',
                                                                            'Cage_No.', 'Animal_No.']




    _metadata.reset_index(inplace=True, drop=True)
    data.reset_index(inplace=True, drop=True)
    features = pd.concat([_metadata, data], axis=1)

    features.index.name = 'scanID'

    print(features)

    print(str(_dir.parent / "full_results.csv"))

    features.to_csv(str(_dir.parent / "full_results.csv"))

    feature_reduction.main(features, org = None, rad_file_path = Path(_dir.parent / "full_results.csv"))



def test_BQ_mach_learn_batch_sp():
    _dir = Path("E:/220913_BQ_tsphere/inputs/features/")

    file_names = [spec for spec in common.get_file_paths(folder=_dir, extension_tuple=".csv")]
    file_names.sort()

    data = [pd.read_csv(spec, index_col=0).dropna(axis=1) for spec in file_names]

    data = pd.concat(
        data,
        ignore_index=False, keys=[os.path.splitext(os.path.basename(spec))[0] for spec in file_names],
        names=['specimen', 'label'])

    data['specimen'] = data.index.get_level_values('specimen')

    _metadata = data['specimen'].str.split('_', expand=True)



    _metadata.columns = ['Date', 'Exp', 'Contour_Method', 'Tumour_Model', 'Position', 'Age',
                                                                            'Cage_No.', 'Animal_No.']




    _metadata.reset_index(inplace=True, drop=True)
    data.reset_index(inplace=True, drop=True)
    features = pd.concat([_metadata, data], axis=1)

    features.index.name = 'scanID'

    print(features)

    print(str(_dir.parent / "full_results.csv"))

    features.to_csv(str(_dir.parent / "full_results.csv"))

    feature_reduction.main(features, org = None, rad_file_path = Path(_dir.parent / "full_results.csv"), batch_test=True)


def test_BQ_concat_batch():
    _dir = Path("Z:/jcsmr/ROLab/Experimental data/Radiomics/Workflow design and trial results/Kyle Drover analysis/220617_BQ_norm_stage_full/sub_normed_features.csv")
    #_dir = Path("E:/220913_BQ_tsphere/inputs/features/")

    # file_names = [spec for spec in common.get_file_paths(folder=_dir, extension_tuple=".csv")]
    # file_names.sort()
    #
    # data = [pd.read_csv(spec, index_col=0).dropna(axis=1) for spec in file_names]
    #
    # data = pd.concat(
    #     data,
    #     ignore_index=False, keys=[os.path.splitext(os.path.basename(spec))[0] for spec in file_names],
    #     names=['specimen', 'label'])
    #
    # data['specimen'] = data.index.get_level_values('specimen')
    #
    # _metadata = data['specimen'].str.split('_', expand=True)
    #
    #
    #
    # _metadata.columns = ['Date', 'Exp', 'Contour_Method', 'Tumour_Model', 'Position', 'Age',
    #                                                                         'Cage_No.', 'Animal_No.']
    #
    #
    #
    #
    # _metadata.reset_index(inplace=True, drop=True)
    # data.reset_index(inplace=True, drop=True)
    # features = pd.concat([_metadata, data], axis=1)
    #
    # features.index.name = 'scanID'
    #
    # print(features)
    #
    # print(str(_dir.parent / "full_results.csv"))
    #
    # features.to_csv(str(_dir.parent / "full_results.csv"))

    features = pd.read_csv(_dir)
    features = features[features.columns.drop(list(features.filter(regex="diagnostics")))]
    features.drop(["scanID"], axis=1, inplace=True)
    feature_reduction.main(features, org = None, rad_file_path = Path(_dir.parent / "full_results.csv"), batch_test=True)






@pytest.mark.skip
def test_feat_reduction():
    feature_reduction.main()

def test_mach_learn_pipeline():
    lama_machine_learning.ml_job_runner("E:/230129_bq_tester/norm_methods/", n_sample=True)


@pytest.mark.skip
def test_radiomic_org_plotting():
    _dir = Path("E:/220607_two_way/g_by_back_data/radiomics_output/features/")

    file_names = [spec for spec in common.get_file_paths(folder=_dir, extension_tuple=".csv")]
    file_names.sort()

    data = [pd.read_csv(spec, index_col=0).dropna(axis=1) for spec in file_names]

    abnormal_embs = ['22300_e8', '22300_e6', '50_e5']

    for i, df in enumerate(data):
        df.index.name = 'org'
        df.name = str(file_names[i]).split(".")[0].split("/")[-1]
        df['genotype'] = 'HET' if 'het' in str(file_names[i]) else 'WT'
        df['background'] = 'C56BL6N' if (('b6ku' in str(file_names[i])) | ('BL6' in str(file_names[i]))) else 'C3HHEH'
        df['HPE'] = 'abnormal' if any(map(str(file_names[i]).__contains__, abnormal_embs)) else 'normal'

    data = pd.concat(
        data,
        ignore_index=False, keys=[os.path.splitext(os.path.basename(spec))[0] for spec in file_names],
        names=['specimen', 'org'])

    line_file = _dir.parent / "full_results.csv"

    data.to_csv(line_file)

    #data_subset = data.select_dtypes(include=np.number)

    for i, org in enumerate(data.index.levels[1]):
        fig, ax = plt.subplots(1, 1, figsize=[56, 60])
        #sns.set(font_scale=0.5)
        o_data = data[np.isin(data.index.get_level_values('org'), org)]

        o_data_subset = o_data.select_dtypes(include=np.number)
        #o_data_subset = o_data_subset.apply(lambda x: (x - x.mean()) / x.std(), axis=0)
        o_data_subset = o_data_subset.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

        tsne = TSNE(perplexity=30,
                    n_components=2,
                    random_state=0,
                    early_exaggeration=250,
                    n_iter=1000,
                    verbose=1)

        tsne_results = tsne.fit_transform(o_data_subset.dropna(axis='columns'))

        o_data['tsne-2d-one'] = tsne_results[:, 0]
        o_data['tsne-2d-two'] = tsne_results[:, 1]
        o_data['org'] = o_data.index.get_level_values('org')
        o_data['specimen'] = o_data.index.get_level_values('specimen')

        o_data['condition'] = o_data['genotype'] + "_" + o_data['background']

        fig, ax = plt.subplots()
        o_data = o_data[o_data['condition'] == 'WT_C3HHEH']
        g = sns.lmplot(
            x="tsne-2d-one", y="tsne-2d-two",
            data=o_data,
            # col_order=['WT_C3HHEH','HET_C3HHEH','WT_C57BL6N','HET_C57BL6N'],
            #col='specimen',
            #col_wrap=5,
            hue="specimen",
            palette='husl',
            fit_reg=False)

        def label_point(x, y, val, ax):
            a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
            for i, point in a.iterrows():
                ax.text(point['x'] + .02, point['y'], str(point['val']), fontsize='xx-small')

        label_point(o_data['tsne-2d-one'], o_data['tsne-2d-two'], o_data['specimen'], plt.gca())



        #g.set(ylim=(np.min(o_data['tsne-2d-two']) - 10, np.max(o_data['tsne-2d-two']) + 10),
         #     xlim=(np.min(o_data['tsne-2d-one']) - 10, np.max(o_data['tsne-2d-one']) + 10))
        plt.savefig("E:/220607_two_way/g_by_back_data/radiomics_output/radiomics_2D_tsne_C3H_wt_" + str(org) + ".png")
        plt.close()


def test_get_rad_data_for_perm():
    _dir = Path("E:/221122_two_way/g_by_back_data/radiomics_output")

    wt_dir = Path("E:/221122_two_way/g_by_back_data/baseline")

    mut_dir = Path("E:/221122_two_way/g_by_back_data/mutants")

    treat_dir = Path("E:/221122_two_way/g_by_back_data/treatment")

    inter_dir = Path("E:/221122_two_way/g_by_back_data/mut_treat")

    results = get_radiomics_data(_dir, wt_dir, mut_dir, treat_dir, inter_dir)

    results.to_csv(str(_dir/"test_dataset.csv"))

def test_permutation_stats():
    """
    Run the whole permutation based stats pipeline.
    Copy the output from a LAMA registrations test run, and increase or decrease the volume of the mutants so we get
    some hits

    """
    lama_permutation_stats.run(stats_cfg)

def test_sns_clustermap():
    url = "https://raw.githubusercontent.com/dorkylever/LAMA/master/lama/tests/clustermap_data.csv"
    X = pd.read_csv(url, index_col=0)
    X.dropna(how='all', inplace=True)
    X.fillna(1, inplace=True)
    # replace missing values with 0
    # replace infinite values with 0

    print("Missing values:", X.isnull().sum().sum())
    print("Infinite values:", np.isposinf(X).sum().sum() + np.isneginf(X).sum().sum())

    # mean_data = X.apply(np.mean, axis=1)

    # std_data = X.apply(np.std, axis=1)

    # constant_rows = X.apply(lambda row: row.nunique() == 1, axis=1)

    # X = X[~constant_rows]

    # na_mean = mean_data[mean_data.isna().any()]

    # na_std = std_data[std_data.iszero().any()]

    cg = sns.clustermap(X,
                        metric="euclidean",
                        cmap=sns.diverging_palette(250, 15, l=70, s=400, sep=1, n=512, center="light",
                                                   as_cmap=True),
                        cbar_kws={'label': 'mean volume ratio'}, square=True,
                        center=1,
                        figsize=[30, len(X) * 0.3])

    # Calculate the mean and standard deviation of each variable
    # X.columns = [col.rsplit('_', 3)[-1] for col in X.columns]

    plt.tight_layout()

    plt.savefig("E:/221122_two_way/permutation_stats/rad_perm_output/two_way/easy_fig.png")
    plt.close()

    X = X.to_numpy()

    print(np.isnan(X).any() | np.isinf(X).any())
    #
    print(X)
    mu = np.mean(X, axis=0)
    sigma = np.std(X, axis=0)
    #
    #
    print(np.all(sigma))
    #
    # # Calculate the z-score matrix
    Z = (X - mu) / sigma
    print(Z)
    #
    # # Calculate the Euclidean distance matrix
    d = np.zeros((Z.shape[0], Z.shape[0]))
    for i in range(Z.shape[0]):
        for j in range(Z.shape[0]):
            d[i, j] = np.sqrt(np.sum((Z[i, :] - Z[j, :]) ** 2))
    #
    print(d)
    print(np.isnan(d).any() | np.isinf(d).any())