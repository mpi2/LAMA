import pandas as pd
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import common

# The maximum allowed size difference (percent) between the extremes of the WT range and the mutant range
MAX_DIFF_MUTANT_SMALLEST_WT = 0.05
MIN_WTS = 8  # The minimum number of baselines needed for statistical analysis


class BaselineSelector(object):
    """
    Given two staging csv files previously created by lama,
            eg:
            -----------
            wt1,700
            wt2,710
            wt3,720
            wt4,730
            ...
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

    def __init__(self, wt_staging: pd.DataFrame, mut_staging: pd.DataFrame, littermate_basenames: list=None, mut_ids: list=None):
        """
        
        Parameters
        ----------
        wt_staging: pd.Dataframe
        mut_staging: pd.Dataframe
        littermate_basenames
        mut_ids: 
            list: if only using a subset of the mutants
            None: if using all the mutants
        """
        if littermate_basenames:
            self.littermate_basenames = [os.path.basename(x) for x in littermate_basenames]
        else:
            self.littermate_basenames = None

        self.wt_df = wt_staging
        self.mut_df = mut_staging
        self.littermate_df = self._get_littermate_dfs()

        # Check staging csv headers
        wt_heads = [x in self.wt_df.columns for x in ['vol', 'value']]
        mut_heads = [x in self.mut_df.columns for x in ['vol', 'value']]

        if not all(wt_heads + mut_heads):
            raise common.LamaDataException("The staging files must contain the headers: vol, value")

        self.wt_df.set_index(self.wt_df['vol'], inplace=True)
        self.wt_df.vol = self.wt_df.vol.astype(str)  # in case the name is purley numeric
        self.mut_df.set_index(self.mut_df['vol'], inplace=True)

        self.mut_ids = mut_ids  # extension-stripped specimen ids
        self.sorted_df = self.wt_df.sort_values(by='value', ascending=True)

        self.excluded_mutants = self._get_mutants_outside_staging_range()

        self.df_filtered_wts = self._generate()

    def _get_littermate_dfs(self):
        if self.littermate_basenames:
            litter_mate_df = self.mut_df[self.mut_df['vol'].isin(self.littermate_basenames)]
            mut_df = self.mut_df[~self.mut_df['vol'].isin(self.littermate_basenames)]
        else:
            litter_mate_df = None
        return litter_mate_df

    def filtered_wt_ids(self, ignore_constraint=False):
        """
        Get the final list of baseline ids to use

        Parameters
        ----------
        ignore_constraint: bool
            if True, do not select by staging criteria, and return all baselines

        Returns
        -------
        list of specimen ids
        """

        if self.df_filtered_wts is not None and not ignore_constraint:
            wt_ids = [str(x) for x in list(self.df_filtered_wts.index)]
            return wt_ids
        elif ignore_constraint:
            return self._generate_without_constraint()
        else:
            return None

    def _get_mutants_outside_staging_range(self):
        """
        Find mutants that are too small (usually) or too large(not seen yet) to analyse as there are no suitably
        stage-matched baselines

        Returns
        -------
        list
            mutant ids to exclude due to extreme of size
        """
        result = []

        wt_min = self.wt_df['value'].min()
        wt_min = wt_min - (wt_min * MAX_DIFF_MUTANT_SMALLEST_WT)
        wt_max = self.wt_df['value'].max()
        wt_max = wt_max + (wt_max * MAX_DIFF_MUTANT_SMALLEST_WT)

        for i, row in self.mut_df.iterrows():
            staging_metric = row['value']

            if wt_min <= staging_metric <= wt_max:
                continue
            else:
                result.append(row.vol)
        return result

    def littermates_to_include(self):
        """
        Get the littermate IDs that are to be included with the baseline set. Littermates that are too large will
        not be included

        Returns
        -------
        list[baseline ids]
        """
        if self.littermate_df is None:
            return None  # There are no littermates to consider

        mut_min = self.mut_df['value'].min()
        mut_max = self.mut_df['value'].max()

        to_include = self.littermate_df[(self.littermate_df['value'] >= mut_min)
                                              & (self.littermate_df['value'] <= mut_max)]
        if len(to_include) == 0:
            return None
        else:
            return common.specimen_ids_from_paths(to_include['vol'].tolist())

    def get_mut_crls(self):
        mut_crls = dict(list(zip(self.mut_df.vol, self.mut_df['value'])))
        return mut_crls

    def get_wt_crls(self):
        wt_crls = dict(list(zip(self.df_filtered_wts.vol, self.df_filtered_wts['value'])))
        return wt_crls

    def _generate_without_constraint(self):
        """
        Generate a list of baseline specimens that are closet to a mutant set in some developmental stage proxy
        Return the nearest regardless of how much they deviate from the range of mutants
        Returns
        -------

        """
        # Assuming we are using this due to the fact that all mutants are larger or smaller than available baselines.
        # So we just pick the min number off the top or bottom of the pile
        return self.sorted_df.vol.tolist()
        wt_min = self.sorted_df.iloc[0].value
        mut_min = self.mut_df['value'].min()
        if mut_min > wt_min:
            result = self.sorted_df[-MIN_WTS:].vol
        else:
            result = self.sorted_df[0: MIN_WTS].vol
        return list(result)

    def _generate(self, max_extra_allowed=MAX_DIFF_MUTANT_SMALLEST_WT):
        """

        Parameters
        ----------

        Returns
        -------

        """
        # If we have a bunch of littermate ids, drop these from the DataFrame used for stage calculation
        if self.littermate_basenames:
            to_drop = []
            for id_ in self.mut_df.vol:
                if common.strip_img_extension(id_) in self.littermate_basenames or id_ in self.littermate_basenames:
                    to_drop.append(id_)
            self.mut_df.drop(to_drop, inplace=True)

        # Only keeps ids specifed in self.mut_ids if self.mut_ids is not None
        # For example there may be hets we don't want to include
        if self.mut_ids:
            for v in self.mut_df.vol:
                if common.strip_img_extension(v) not in self.mut_ids:
                    self.mut_df = self.mut_df[self.mut_df.vol != v]

        # remove excluded mutants
        for v in self.mut_df.vol:
            if common.strip_img_extension(v) in self.excluded_mutants:
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
                if current_min_idx < 0:
                    min_reached = True

                current_max_idx += 1
                if current_max_idx > len(self.sorted_df) - 1:
                    max_reached = True

                # Get the next biggest vol
                if not max_reached:
                    try:
                        nbv = self.sorted_df.iloc[current_max_idx]
                    except IndexError:
                        max_reached = True
                    else:
                        if nbv.value <= expanded_max:
                            new_additions_inices.append(nbv)
                            vol_num += 1
                            if vol_num >= MIN_WTS:
                                break
                        else:
                            max_reached

                # Get the next smallest vol
                if not min_reached:
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



