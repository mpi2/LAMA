"""
Henrik wants plots of the organ volume results hits.


25/02/19
Extending it to do all lines. Move into the lama pipeline at some point

27/02/19
Added scatter plots

"""

from pathlib import Path
import math
from typing import List

import numpy as np
import seaborn as sns
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import pandas as pd
from logzero import logger as logging

from lama.common import getfile_startswith_endswith
from lama.qc import formatting

organ_vol = 'organ volume'   # Y label
wev = 'whole embryo volume'  # x label for scatter plots


def pvalue_dist_plots(null: pd.DataFrame, alt: pd.DataFrame, thresholds: pd.DataFrame, outdir: Path):
    """
    Generate a series of histograms containing null and alternative distribution overlaid.
    Create a vertical line where the p-value threshold was set

    Parameters
    ----------
    null, alt
        each column shoul dbe a label number. Each row a p-value from a single permutation test
    thresholds
        index: label (number)
        must also have column 'p_thresh'
    outdir
        where to put the plots
    """

    def hist(values):
        hist, bins = np.histogram(values, 100)
        hist = hist / np.sum(hist)
        width = 1.0 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width, alpha=0.5)

    x_label = 'log(p)'
    alt = np.log(alt)
    null = np.log(null)

    for col in alt:
        try:
            thresh = thresholds.loc[int(col), 'p_thresh']
            log_thresh = np.log(thresh)

            hist(alt[col])
            hist(null[col])
            plt.xlabel(x_label)

            outpath = outdir / f'{col}.png'

            plt.axvline(log_thresh, 0, 1, alpha=0.4, color='g')

            plt.legend(labels=['p threshold = {}'.format(format(thresh, '.3g')), 'alt', 'null'])
            plt.ylabel('Density')
            plt.title(col)
            plt.savefig(outpath)
            plt.close()
        except ValueError:
            logging.warn(f'Skipping pvalue dist plot for {col}')
            continue


def make_plots(mut_lines_dir: Path,
               wt_organ_vols: pd.DataFrame,
               wt_staging: pd.DataFrame,
               label_meta_file: Path,
               stats_root_dir: Path,
               skip_no_analysis= False,
               organ_subset: List = [],
               extra_dir: Path = Path('')):
    """

    Parameters
    ----------
    mut_lines_dir
        Lama registration root. eg: mutants/output  with each sibdir containing a line
    wt_organ_vols
        Aggregated organ volumes for each baseline
    wt_staging
        Aggregated staging info for each baseline
    label_meta_file
        CSV of atlas metadata
    stats_root_dir
        Contains a folder for each line
    skip_no_analysis
        If there is a no_analysis=True flag in the meta data csv, skip this organ
    organ_subset
        plot only the labels with these label numbers. Or can be used to order the output of the plots
    extra_dir
        Bit of a bodge, but if doing the non-permutation-based stats, the organ vol csv is in a directory below.
        Give the name here (currently 'organ_volumes')
    """

    label_meta = pd.read_csv(label_meta_file, index_col=0)

    wt_staging.rename(columns={'line': 'genotype'}, inplace=True)

    for mut_line_dir in mut_lines_dir.iterdir():

        if not mut_line_dir.is_dir():
            continue

        print(mut_line_dir.name)

        stats_line_dir = stats_root_dir / mut_line_dir.name / extra_dir  # extra_dir does nothing if == ''

        line = mut_line_dir.name

        #TODO: Get file by startswith line name and endswith extension (Could be date of analysis in middle)
        # Rather tan just getting any CSVs in there

        stats_result_file = getfile_startswith_endswith(stats_line_dir, line, '.csv')

        # Get mutant staging and organ volumes
        line_vols = []
        line_stage = []

        for spec_dir in mut_line_dir.iterdir():
            if str(spec_dir).endswith('_'):
                continue

            staging = pd.read_csv(spec_dir / 'output' / 'staging_info_volume.csv', index_col=0)
            organ_vols = pd.read_csv(spec_dir / 'output' / 'organ_volumes.csv', index_col=0)

            line_vols.append(organ_vols)
            line_stage.append(staging)

        df_stage_mut = pd.concat(line_stage, axis=0)
        df_stage_mut['genotype'] = 'mutant'
        df_stage_mut.rename(columns={'value': 'staging'},  inplace=True) # Get rid of this
        df_vol_mut = pd.concat(line_vols, axis=0)
        df_hits = pd.read_csv(stats_result_file, index_col=0)

        staging_df = pd.concat([wt_staging, df_stage_mut])
        staging_df.rename(columns={'staging': wev}, inplace=True)
        vol_df = pd.concat([wt_organ_vols, df_vol_mut])




        if 'significant_cal_p' in df_hits:  # 'permutation stats
            hits: pd.DataFrame = df_hits[df_hits['significant_cal_p'] == True]
        elif 'significant_bh_q_5' in df_hits:
            hits: pd.DataFrame = df_hits[df_hits['significant_bh_q_5'] == True]
        else:
            logging.error("Plots not made: Stats output file must have 'significant_cal_p' or 'significant_bh_q_5' column")

        if ('organ_system_name' in label_meta.columns) and ('organ_system_name' not in hits):
            # Sort by organ system if present in atlas metadata
            hits = hits.merge(label_meta[['organ_system_name']], how='left', left_index=True, right_index=True)
            hits.sort_values(by='organ_system_name', inplace=True)

        if skip_no_analysis:
            # Skip organ that are flagged with no_analysis in the atlas metadata file
            if 'no_analysis' not in hits:
                hits = hits[hits['no_analysis'] != True]

        if len(hits) < 1:
            logging.info(f'No hits, so Skipping organ vol plots for: {mut_line_dir.name}')
            continue

        st = wt_staging['staging']
        normed_wt = wt_organ_vols.div(st, axis=0)

        normed_mut = df_vol_mut.div(df_stage_mut['staging'], axis=0)

        numcol = 6 if len(hits) > 5 else len(hits)
        numrows = math.ceil(len(hits) / numcol)

        figsize_y = 5 * numrows
        figsize_x = 5 * numcol

        fig = Figure(figsize=(figsize_x, figsize_y))
        FigureCanvas(fig)

        fig_scat = Figure(figsize=(figsize_x, figsize_y))
        FigureCanvas(fig_scat)

        if organ_subset:
            labels_to_plot = organ_subset

            if len(set(labels_to_plot).intersection(hits.index)) != len(labels_to_plot):
                print('Some label numbers in organ_subset are not in the hits DataFrame')
                return
        else:
            labels_to_plot = hits.index

        # for i, (label, row) in enumerate(hits.iterrows()):
        for i, label in enumerate(labels_to_plot):
            label_name: str = hits.loc[label, 'label_name']
            axes = fig.add_subplot(numrows, numcol, i + 1)
            axes.tick_params(labelsize=18)
            axes.set_yticklabels([])

            label = str(label)

            wt = normed_wt[[label]]
            wt['genotype'] = 'baseline'

            mut = normed_mut[[label]]
            mut['genotype'] = line

            df = pd.concat([wt, mut])
            df.rename(columns={label: organ_vol}, inplace=True)

            min_ = df[organ_vol].min() - (df[organ_vol].min() * 0.1)
            max_ = df[organ_vol].max() + (df[organ_vol].max() * 0.1)



            # sns.boxplot(x="genotype", y="organ volume", data=df, orient='v',
            #             ax=axes, boxprops=boxprops)

            axes.tick_params(labelsize=18)

            ax = sns.swarmplot(x="genotype", y='organ volume', data=df, orient='v',
                             ax=axes)

            ax.set_ylim(min_, max_)

            for patch in ax.artists:
                r, g, b, a = patch.get_facecolor()
                patch.set_facecolor((r, g, b, 0.0))

            if 'short_name' in label_meta:
                label_name = label_meta.at[int(label), 'short_name']
            title = label_name.replace('_', ' ')
            title = title.capitalize()
            ax.set_ylabel('')
            ax.set_xlabel('')

            ax.set_title(title, fontsize=20)

            ###Scatter
            s_axes = fig_scat.add_subplot(numrows, numcol, i + 1)
            s_axes.tick_params(labelsize=18)

            # Get rid of hard-coding
            voxel_size = 27.0
            um3_conv_factor = voxel_size ** 3  # To convert voxels to um3
            um3_to_mm3_conv_factor = 1e9

            scattter_df = staging_df.join(vol_df[[label]]).rename(columns={label: organ_vol})
            scattter_df[organ_vol] = (scattter_df[organ_vol] * um3_conv_factor) / um3_to_mm3_conv_factor
            scattter_df[wev] = (scattter_df[wev] * um3_conv_factor) / um3_to_mm3_conv_factor

            scattter_df['normalised_organ_vol'] = scattter_df[organ_vol] / scattter_df[wev]

            sax = sns.scatterplot(y=organ_vol, x=wev, ax=s_axes, hue='genotype',
                                  data=scattter_df)

            sax.set(xlabel='Whole embryo volume (mm^3)')
            sax.set(ylabel='Organ volume (mm^3)')

            sax.set_title(title, fontsize=16)
            sax.ticklabel_format(style='sci',scilimits=(0, 0))

            # x 10^7 instead of 1e7
            sax.xaxis.major.formatter._useMathText = True
            sax.yaxis.major.formatter._useMathText = True

            formatting.label_offset(sax)

        fig.subplots_adjust(top=0.8)  # TODO fix this for larger plot
        fig.suptitle(line, fontsize=30,  y=0.98)
        # fig.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        if skip_no_analysis:
            box_name = f'{line}_boxplots_no_analysis.png'
        else:
            box_name = f'{line}_boxplots.png'

        # TODO: Fix the boxplot or swarm plot output
        # fig.savefig(stats_line_dir / box_name)

        fig_scat.subplots_adjust(top=0.8, wspace=0.35, hspace=0.4)  # TODO fix this for larger plot
        fig_scat.suptitle(line, fontsize=30,  y=0.98)
        fig_scat.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        if skip_no_analysis:
            scatter_name = f'{line}_scatter_plots_no_analysis_normalised.png'
        else:
            scatter_name = f'{line}_scatter_plots.png'
        fig_scat.savefig(stats_line_dir / scatter_name)





