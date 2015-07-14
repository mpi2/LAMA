#!/usr/bin/python

import csv
from collections import defaultdict
from scipy import stats

def organvolume_stats(wt_csv, mut_csv):
    """"
    Given two csv files with organ volume calculation, calculate some statistics and write report

    Parameters
    ----------
    wt_volumes: str
        path to wildtype csv file
    mut_volumes: str
        path to mut csv files
    """

    wt_data = defaultdict(list)
    mut_data = defaultdict(list)

    _get_label_vols(wt_csv, wt_data)
    _get_label_vols(mut_csv, mut_data)

    _tstats(wt_data, mut_data)


def _tstats(wt_data, mut_data):
    for label, wt_values in wt_data.iteritems():
        mut_values = mut_data[label]

        stats.ttest_ind(mut_values, wt_values)





def _get_label_vols(file_path, output_dict):


    with open(file_path, 'r') as csvfile:
        reader = csv.reader(csvfile)
        header = reader.next()
        for row in reader:
            for label_name, vol_calculation in zip(header[1:], row[1:], ):  # Skip the vol name column
                output_dict[label_name].append(vol_calculation)


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser("Calculate organ volume statistics")
    parser.add_argument('-w', dest='wt', help='wt volume calculations csv', required=True)
    parser.add_argument('-m', dest='mut', help='mut volume calculations csv', required=True)
    args = parser.parse_args()

    organvolume_stats(args.wt, args.mut)