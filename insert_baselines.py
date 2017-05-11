"""
insert_baselines.py

The following scenario was in mind when making this module.

You have a set of baseline data theat has been registered towards your current population average.
It contains, amongst other things:

    registered images
    jobobians
    and poibly inverted labels - which might be needed for calulating staging or organ volumes

This module will insert the above data into the correct locations of the WT test set heirachy
If you have inverted labels that are used to calulate staging metrics, these will need to be run manually after insertion
"""


def insert(lama_output_folder, wt_output_folder, vol_ids, lama_config):
    """
    
    Parameters
    ----------
    lama_results_folder: str
        path to the lama run, where there are baselines located. Normally called 'output'
    wt_folder: str
        output directory of WT test set run - namd 'output'
    vol_ids: list or str
        list of ids (file name basenames)
        str path to csv containing the volume basenames
    lama_config: str
        main lama config file used to generate the reults. This is needed to find the correct paths

    Returns
    -------

    """


if __name__ == '__main__':
    lama_output_folder = '/home/neil/sig/LAMA_results/E14.5/compare_cbx2_male_female_variances/male_female_and_cbx2/all_inputs/mutant_runs/nras/output'
    wt_output_folder = '/home/neil/sig/LAMA_results/E14.5/compare_cbx2_male_female_variances/male_female_and_cbx2/all_inputs/output'
    vol_ids = '/home/neil/sig/LAMA_results/E14.5/compare_cbx2_male_female_variances/male_female_and_cbx2/all_inputs/mutant_runs/nras/liitermates.csv'
    insert()