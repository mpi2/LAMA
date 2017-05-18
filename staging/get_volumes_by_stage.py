import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import common


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
        csv path with staging info (sacling factors for each id for example)
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

    def __init__(self, wt_staging_file, mut_staging_file, littermate_basenames=None):
        self.littermate_basenames = littermate_basenames
        self.wt_df = pd.read_csv(wt_staging_file)
        self.wt_df.set_index(self.wt_df['vol'], inplace=True)
        self.mut_df = pd.read_csv(mut_staging_file)
        self.mut_df.set_index(self.mut_df['vol'], inplace=True)
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

    def filtered_wt_ids(self):
        if self.df_filtered_wts is not None:
            return list(self.df_filtered_wts.index)
        else:
            return None

    def _generate(self, max_extra_allowed=0.08):
        """

        Parameters
        ----------

        Returns
        -------

        """
        if self.littermate_basenames:
            self.mut_df.drop(self.littermate_basenames, inplace=True)
        min_wts = 8

        mut_min = self.mut_df['value'].min()
        mut_max = self.mut_df['value'].max()

        sorted_df = self.wt_df.sort(columns=['value'], ascending=True)

        # First off, get all WT volumes within the range of the mutants
        filtered_df = sorted_df[sorted_df['value'].between(mut_min, mut_max, inclusive=True)]

        if len(filtered_df) < min_wts:
            # Now try expanding the allowed range of WTs to see if we then have enough
            expanded_min = mut_min - (mut_min * max_extra_allowed)
            expanded_max = mut_max + (mut_min * max_extra_allowed)
            filtered_df = sorted_df[sorted_df['value'].between(expanded_min, expanded_max, inclusive=True)]

            if len(filtered_df) >= min_wts:
                # trim off either side of the range until we get 8. Too many WTs outside the muta range may skey results
                while len(filtered_df) > min_wts:
                    filtered_df.drop(filtered_df[0], inplace=True)  # take one of the start
                    if len(filtered_df) <= min_wts:
                        break
                    filtered_df.drop(filtered_df[:-1], inplace=True)  # take one of the end
                    if len(filtered_df) <= min_wts:
                        break
            else:
                return None
        # Return the staged list of ids
        return filtered_df



