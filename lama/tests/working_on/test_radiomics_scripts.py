
from pathlib import Path
from lama.radiomics.radiomics import radiomics_job_runner
from lama import common
import os
from lama.common import cfg_load
from lama.img_processing import normalise
import logging
import seaborn as sns
import pandas as pd
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pytest
import numpy as np

@pytest.mark.skip
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
        df['background'] = 'C56BL6N' if (('b6ku' in str(file_names[i]))|('BL6' in str(file_names[i]))) else 'C3HHEH'

        df['HPE'] = 'abnormal' if any(map(str(file_names[i]).__contains__, abnormal_embs)) else 'normal'

    data = pd.concat(
        data,
        ignore_index=False, keys=[os.path.splitext(os.path.basename(spec))[0] for spec in file_names],
        names=['specimen', 'org'])

    line_file = _dir.parent / "full_results.csv"

    data.to_csv(line_file)

    data_subset = data.select_dtypes(include=np.number)

    data_subset = data_subset.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    data_subset = data_subset.apply(lambda x: (x - x.mean()) / x.std(), axis=0)


    tsne = TSNE(perplexity=30,
                n_components=2,
                random_state=0,
                early_exaggeration=250,
                n_iter=1000,
                verbose=1)

    #print(data_subset.dropna(axis='columns'))

    tsne_results = tsne.fit_transform(data_subset.dropna(axis='columns'))

    color_class = data.index.get_level_values('org')

    # fig, ax = plt.subplots(figsize=[55, 60])
    # cluster.tsneplot(score=tsne_results, show=True, theme='dark', colorlist=color_class)

    data['tsne-2d-one'] = tsne_results[:, 0]
    data['tsne-2d-two'] = tsne_results[:, 1]
    data['org'] = data.index.get_level_values('org')
    data['specimen'] = data.index.get_level_values('specimen')
    data['condition'] = data['genotype'] + "_" + data['background']



    # data['tsne-4d-three'] = tsne_results[:, 1]

    # fig = plt.figure()
    # ax = fig.add_subplot(110, projection='3d')

    # plt.figure(figsize=(29, 30))
    # data["org"] = data["org"].astype(str)
    # cmap = ListedColormap(sns.color_palette("husl", 255).as_hex())
    # colours = data.index.get_level_values('org')

    # ax.scatter2D(data['tsne-3d-one'], data['tsne-3d-two'], data['tsne-3d-three'], c = colours, cmap=cmap)
    # plt.show(

    fig, ax = plt.subplots(figsize=[56, 60])
    data = data[data['condition'] == 'WT_C3HHEH']
    g = sns.lmplot(
        x="tsne-2d-one", y="tsne-2d-two",
        data=data,
        #col_order=['WT_C3HHEH','HET_C3HHEH','WT_C57BL6N','HET_C57BL6N'],
        col='org',
        col_wrap=5,
        hue="org",
        palette='husl',
        fit_reg=False)


    g.set(ylim=(np.min(data['tsne-2d-two'])-10, np.max(data['tsne-2d-two'])+10),
          xlim=(np.min(data['tsne-2d-one'])-10, np.max(data['tsne-2d-one'])+10))
    plt.savefig("E:/220607_two_way/g_by_back_data/radiomics_output/features/radiomics_2D_tsne_C3H_wt_v2_org.png")

    #plt.close()

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


