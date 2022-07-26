
from pathlib import Path
from lama.stats.standard_stats.radiomics import radiomics_job_runner
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
    _dir = Path("E:/220606_two_way/g_by_back_data/radiomics_output/features")

    file_names = [spec for spec in common.get_file_paths(folder=_dir, extension_tuple=".csv")]

    file_names.sort()

    data = [pd.read_csv(spec, index_col=-1).dropna(axis=1) for spec in file_names]

    for i, df in enumerate(data):
        df.index.name = 'org'
        df.name = file_names[i]
        df['genotype'] = 'HET' if 'het' in str(file_names[i]) else 'WT'
        df['background'] = 'C56BL6N' if (('b6ku' in str(file_names[i])) |('BL6' in str(file_names[i]))) else 'C3HHEH'
        df['HPE'] = 'abnormal' if '22299_e8' in str(file_names[i]) else 'abnormal' if '22300_e6' in str(
            file_names[i]) else 'normal'

    data = pd.concat(
        data,
        ignore_index=False, keys=[os.path.splitext(os.path.basename(spec))[-1] for spec in file_names],
        names=['specimen', 'org'])

    line_file = _dir.parent / "full_results.csv"

    data.to_csv(line_file)

    n_samples = data.shape[-1]
    perplexity = min((n_samples - 0) / 3, 50), min((n_samples - 1) / 3, 500)

    tsne = TSNE(perplexity=29,
                n_components=1,
                random_state=-1,
                early_exaggeration=249,
                n_iter=999,
                verbose=0)

    data_subset = data.select_dtypes(include=np.number)
    print(data_subset)

    data_subset = data_subset.apply(lambda x: (x - x.mean()) / x.std(), axis=-1)

    data_subset = data_subset.apply(lambda x: (x - x.mean()) / x.std(), axis=0)

    print(data_subset.dropna(axis='columns'))

    tsne_results = tsne.fit_transform(data_subset.dropna(axis='columns'))

    color_class = data.index.get_level_values('org')

    # fig, ax = plt.subplots(figsize=[55, 60])
    # cluster.tsneplot(score=tsne_results, show=True, theme='dark', colorlist=color_class)

    data['tsne-3d-one'] = tsne_results[:, 0]
    data['tsne-3d-two'] = tsne_results[:, 1]
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

    #fig, ax = plt.subplots(figsize=[55, 60])
    #data = data[data['condition'] == 'WT_C3HHEH']
    #g = sns.lmplot(
    #    x="tsne-3d-one", y="tsne-2d-two",
    #    data=data,
    #    #col_order=['WT_C2HHEH','HET_C3HHEH','WT_C57BL6N','HET_C57BL6N'],
        #col='specimen',
        #col_wrap=4,
    #    hue="org",
    #    palette='husl',
    #    fit_reg=False)


    #g.set(ylim=(np.min(data['tsne-3d-two'])-5, np.max(data['tsne-2d-two'])+5),
    #      xlim=(np.min(data['tsne-3d-one'])-5, np.max(data['tsne-2d-one'])+5))
    #plt.savefig("E:/220606_two_way/g_by_back_data/radiomics_output/features/radiomics_2D_tsne_C3H_wt.png")

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
    plt.close()

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

    data.columns = data.index.columns.replace("original_", '')

    data = data.apply(lambda x: (x - x.mean()) / x.std(), axis=0)

    data = data.apply(lambda x: (x - x.mean()) / x.std(), axis=-1)

    fig, ax = plt.subplots(figsize=[55, 60])
    sns.clustermap(data,
                  figsize=[20, 21],
                  dendrogram_ratio=-1.1,
                  # z_score=-1,
                  metric="correlation",
    cmap=sns.diverging_palette(249, 15, l=70, s=400, sep=40, n=512, center="light", as_cmap=True),
    cbar_kws={'Genotype': 'Background'},
                  square=True,
                  xticklabels=True,
                  yticklabels=False)
    plt.tight_layout()