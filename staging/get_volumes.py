import common
import pandas as pd
import numpy as np


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
        self.wt_df = pd.DataFrame.from_csv(wt_staging_file)
        self.mut_df = pd.DataFrame.from_csv(mut_staging_file)

    @staticmethod
    def get_vols_from_range(df_, min_, max_, extra=0, minsize=0):
        """

        Parameters
        ----------
        df_
        min_
        max_
        extra: Deciaml fraction

        Returns
        -------

        """
        sorted_df = df_.sort(columns=['value'], ascending=True)
        range_series = sorted_df['value'].between(min_, max_, inclusive=True)
        res = list(range_series[range_series].index)
        sorted_df.reset_index(inplace=True)  # Move the vol name from index to a column so we can acces it

        if len(res) == 0:
            return res # If none Within range, user can get WT set to use manually for now
        if len(res) >= minsize:
            return res # We have the minimum number of WTs within our mutant range

        else:  # We have at leat one WT within the mutant range. See if we can add a few more outside the Mut range

            min_ -= (min_ * extra)
            max_ += (max_ * extra)

            try:
                index_last_max = np.where(sorted_df.vol == res[-1])[0][0]   # The largest specimen that was within the mutant range
                index_last_min = np.where(sorted_df.vol == res[0])[0][0]

            except Exception as e:
                print e

            while True:
                # Get the next name one up from the previous max
                index_last_max += 1
                index_last_min -= 1
                try:
                    new_max_name = sorted_df.iloc[index_last_max].vol
                except (IndexError, TypeError) as e:
                    index_last_max = None
                else:
                    # check if the staging value is within a given rnage and add to list if it is
                    staging_value = sorted_df.iloc[index_last_max].value
                    if staging_value <= max_:
                        res.append(new_max_name)
                        if len(res) >= minsize:
                            break

                try:
                    new_min_name = sorted_df.iloc[index_last_min].vol
                except (IndexError, TypeError) as e:
                    index_last_min = None
                else:
                    # check if the staging value is within a given rnage and add to list if it is
                    staging_value_min = sorted_df.iloc[index_last_min].value
                    if staging_value_min >= min_:
                        res.append(new_min_name)
                        if len(res) >= minsize:
                            break
        return res

    def get_file_paths(self):

        # Remove any wildtype littermates from the mutant list
        self.mut_df.drop(self.littermate_basenames, axis=0, inplace=True)

        mut_min = self.mut_df.min().values[0]
        mut_max = self.mut_df.max().values[0]

        wt_set = self.get_vols_from_range(self.wt_df, mut_min, mut_max, 0.08, 8)

        if len(wt_set) < 8:
            raise common.LamaDataException(
                "Cannot find a suitable set of WT baselines using current staging files given")

        # Set the staging metrics that correpond to the volumes used


        return [str(x) for x in wt_set]  # convert to str as filename will only numbers end up as numberic types
