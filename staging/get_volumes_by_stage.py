import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
from os.path import splitext
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import common

# The maximum allowed size difference between the extremes of the WT range and the mutant range
MAX_PERCENT_LARGER = 0.15
MIN_WTS = 8


class VolumeGetter(object):
    """
    Given two staging csv files previously created by lama,
            eg:
            -----------
            wt1.nrrd,700
            wt2.nrrd,710
            wt3.nrrd,720
            wt4.nrrd,730....
            ------------
    get the list of wts that are nearest to the range of the mutants
    Parameters
    ----------
    wt_staging_file: str
        csv path with staging info (Affine scaling factors for each id for example)
    mut_staging_file
        csv path with staging info
    littermate_basenames: list
        Staging should exclude the size range of any wild type littermates as they are quite often a bit larger
        This is a list of those to exclude from the staging calculation

    Returns
    -------
    list of wild type specimen ids to use
    None if no suitable range of baselines could be found

    """

    def __init__(self, wt_staging_file, mut_staging_file, littermate_basenames=None, mut_ids=None):
        """
        
        Parameters
        ----------
        wt_staging_file
        mut_staging_file
        littermate_basenames
        mut_ids: 
            list: if only using a subset of the mutants
            None: if using all the mutants
        """
        self.littermate_basenames = littermate_basenames
        self.wt_df = pd.read_csv(wt_staging_file)
        self.wt_df.set_index(self.wt_df['vol'], inplace=True)
        self.mut_df = pd.read_csv(mut_staging_file)
        self.mut_df.set_index(self.mut_df['vol'], inplace=True)
        self.mut_ids = mut_ids  # extension-stripped specimen ids
        self.sorted_df = self.wt_df.sort_values(by='value', ascending=True)

        self.df_filtered_wts = self._generate()

    def plot(self, wt_label='wt', mut_label='mutant', outpath=None):

        x = [wt_label] * len(self.df_filtered_wts)
        x += [mut_label] * len(self.mut_df)
        y = list(self.df_filtered_wts['value']) + list(self.mut_df['value'])
        df = pd.DataFrame.from_dict({'genotype': x, 'size': y}, orient='index').T
        df['size'] = df['size'].astype('float')
        ax = plt.axes()
        sns.swarmplot(x='genotype', y='size', data=df)
        ax.set_title("Staging metrics")
        if outpath:
            plt.savefig(outpath)
        else:
            plt.show()

    def filtered_wt_ids(self, ignore_constraint=False):
        if self.df_filtered_wts is not None and not ignore_constraint:
            return [str(x) for x in list(self.df_filtered_wts.index)]
        elif ignore_constraint:
            return self._generate_without_constraint()
        else:
            return None

    def _generate_without_constraint(self):
        """
        Generate a list of baseline specimens that are closet to a mutant set in some developmental stage proxy
        Return the nearest regardless of how much they deviate from the range of mutants
        Returns
        -------

        """
        # Assuming we are using this due to the fact that all mutants are larger or smaller than available baselines.
        # So we just pick the min number off the top or bottom of the pile
        wt_min = self.sorted_df.iloc[0].value
        mut_min = self.mut_df['value'].min()
        if mut_min > wt_min:
            result = self.sorted_df[-MIN_WTS:].vol
        else:
            result = self.sorted_df[0: MIN_WTS].vol
        return list(result)

    def _generate(self, max_extra_allowed=MAX_PERCENT_LARGER):
        """

        Parameters
        ----------

        Returns
        -------

        """
        if self.littermate_basenames:
            to_drop = []
            for id_ in self.mut_df.vol:
                if common.strip_extension(id_) in self.littermate_basenames or id_ in self.littermate_basenames:
                    to_drop.append(id_)
            self.mut_df.drop(to_drop, inplace=True)

        if self.mut_ids:
            for v in self.mut_df.vol:
                if common.strip_extensions([v])[0] not in self.mut_ids:
                    self.mut_df = self.mut_df[self.mut_df.vol != v]

        mut_min = self.mut_df['value'].min()
        mut_max = self.mut_df['value'].max()

        # First off, get all WT volumes within the range of the mutants
        filtered_df = self.sorted_df[self.sorted_df['value'].between(mut_min, mut_max, inclusive=True)]

        if len(filtered_df) < MIN_WTS:
            if len(filtered_df) < 1:
                # No Wildtypes withn the range of the mutants. So the mutants must all be larger or smaller
                # return None for now
                return None
            else:
                # This is a bodge until I can understand Pandas better
                volnames = list(self.sorted_df.vol)
                min_vol_name = filtered_df.iloc[0].vol
                max_vol_name = filtered_df.iloc[-1].vol
                current_min_idx = volnames.index(min_vol_name)
                current_max_idx = volnames.index(max_vol_name)
            # Now try expanding the allowed range of WTs to see if we then have enough
            expanded_min = mut_min - (mut_min * max_extra_allowed)
            expanded_max = mut_max + (mut_min * max_extra_allowed)
            new_additions_inices = []
            vol_num = filtered_df.vol.size
            min_reached = False
            max_reached = False

            while vol_num < MIN_WTS:

                current_min_idx -= 1
                current_max_idx += 1
                # Get the next biggest vol
                try:
                    nbv = self.sorted_df.iloc[current_max_idx]
                except IndexError:
                    max_reached = True
                if nbv.value <= expanded_max:
                    new_additions_inices.append(nbv)
                    vol_num += 1
                    if vol_num >= MIN_WTS:
                        break
                else:
                    max_reached

                # Get the next smallest vol
                try:
                    nsv = self.sorted_df.iloc[current_min_idx]
                except IndexError:
                    min_reached = True
                if nsv.value >= expanded_min:
                    new_additions_inices.append(nsv)
                    vol_num += 1
                    if vol_num >= MIN_WTS:
                        break
                else:
                    min_reached

                if all([max_reached, min_reached]):
                    return None
        else:
            return filtered_df
        # Return the staged list of ids
        for row in new_additions_inices:
            filtered_df = filtered_df.append(row)
        return filtered_df



