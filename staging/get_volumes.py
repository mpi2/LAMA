import common
import pandas as pd


class VolumeGetter(object):
    """

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
        range_series = df_['value'].between(min_, max_)
        res = list(range_series[range_series].index)
        if len(res) == 0:
            return res
        if res > minsize:
            return res

        else:

            min_ -= (min_ * extra)
            max_ += (max_ * extra)



    def get_file_paths(self):
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
        # Remove any wildtype littermates from the mutant list
        self.mut_df.drop(self.littermate_basenames, inplace=True)

        mut_min = self.mut_df.min().values[0]
        mut_max = self.mut_df.max().values[0]

        wt_set = self.get_vols_from_range(self.wt_df, mut_min, mut_max, 0.08, 8)

        if len(wt_set) < 8:
            raise common.LamaDataException(
                "Cannot find a suitable set of WT baselines using current staging files given")

        return [str(x) for x in wt_set]  # convert to str as filename will only numbers end up as numberic types
