"""
Henrik wants plots of the organ volume results hits.


25/02/19
Extending it to do all lines. Move into the lama pipeline at some point

27/02/19
Added scatter plots

"""

from pathlib import Path
import math

import seaborn as sns
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import pandas as pd

from lama.common import getfile_endswith


def boxplotter(mut_lines_dir: Path, wt_organ_vols, wt_staging, label_meta_file, stats_dir):

    label_meta = pd.read_csv(label_meta_file, index_col=0)

    for mut_line_dir in mut_lines_dir.iterdir():

        if not mut_line_dir.is_dir():
            continue

        print(mut_line_dir.name)
        stats_line_dir = stats_dir / mut_line_dir.name

        line = mut_line_dir.name

        stats_result_file = getfile_endswith(stats_line_dir, '.csv')

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
        vol_df = pd.concat([wt_organ_vols, df_vol_mut])

        hits: pd.DataFrame = df_hits[df_hits['significant_cal_p'] == True]

        if len(hits) < 1:
            continue

        hits = hits.merge(label_meta[['organ_system_name']], how='inner', left_index=True, right_index=True)
        hits.sort_values(by='organ_system_name', inplace=True)

        st = wt_staging['staging']
        normed_wt = wt_organ_vols.div(st, axis=0)

        normed_mut = df_vol_mut.div(df_stage_mut['staging'], axis=0)
        #
        # if line != '1200014J11RIK':
        #     return

        numcol = 6 if len(hits) > 5 else len(hits)
        numrows = math.ceil(len(hits) / numcol)

        figsize_y = 5 * numrows
        figsize_x = 5 * numcol

        fig = Figure(figsize=(figsize_x, figsize_y))
        FigureCanvas(fig)

        fig_scat = Figure(figsize=(figsize_x, figsize_y))
        FigureCanvas(fig_scat)

        i = 0  # controls the axes placement

        boxprops = dict(linewidth=1, facecolor=None)

        for i, (label, row) in enumerate(hits.iterrows()):

            axes = fig.add_subplot(numrows, numcol, i + 1)
            axes.tick_params(labelsize=18)
            axes.set_yticklabels([])

            label_name: str = row['label_name']

            label = str(label)

            wt = normed_wt[[label]]
            wt['genotype'] = 'baseline'

            mut = normed_mut[[label]]
            mut['genotype'] = line

            df = pd.concat([wt, mut])
            df.rename(columns={label: label_name}, inplace=True)

            min_ = df[label_name].min() - (df[label_name].min() * 0.1)
            max_ = df[label_name].max() + (df[label_name].max() * 0.1)

            sns.boxplot(x="genotype", y=label_name, data=df, orient='v',
                        ax=axes, boxprops=boxprops)

            axes.tick_params(labelsize=18)

            ax = sns.swarmplot(x="genotype", y=label_name, data=df, orient='v',
                             ax=axes)

            ax.set_ylim(min_, max_)

            for patch in ax.artists:
                r, g, b, a = patch.get_facecolor()
                patch.set_facecolor((r, g, b, 0.0))

            title = label_name.replace('_', ' ')
            title = title.capitalize()
            ax.set_ylabel('')
            ax.set_xlabel('')

            ax.set_title(title, fontsize=20)

            ###~~Scatter
            s_axes = fig_scat.add_subplot(numrows, numcol, i + 1)
            s_axes.tick_params(labelsize=18)
            s_axes.set_yticklabels([])
            # s_axes.set_p

            sax = sns.scatterplot(y=vol_df[label], x=staging_df['staging'], ax=s_axes, hue=df['genotype'])
            sax.set_title(title, fontsize=20)

        fig.subplots_adjust(top=0.8)  # TODO fix this for larger plot
        fig.suptitle(line, fontsize=30,  y=0.98)
        fig.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        fig.savefig(stats_line_dir / f'{line}_boxplots.png')

        fig_scat.subplots_adjust(top=0.8)  # TODO fix this for larger plot
        fig_scat.suptitle(line, fontsize=30,  y=0.98)
        fig_scat.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        fig_scat.subplots_adjust(hspace=0.4)
        fig_scat.savefig(stats_line_dir / f'{line}_scatter_plots.png')





