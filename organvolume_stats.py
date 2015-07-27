#!/usr/bin/python

import csv
from collections import defaultdict

from scipy import stats
import numpy as np

def organvolume_stats(wt_csv, mut_csv, outfile):
    """"
    Given two csv files with organ volume calculation, calculate some statistics and write report

    Parameters
    ----------

        path to wildtype csv file
    mut_volumes: str
        path to mut csv files
    """

    wt_data = defaultdict(list)
    mut_data = defaultdict(list)

    _get_label_vols(wt_csv, wt_data)
    _get_label_vols(mut_csv, mut_data)

    results = _tstats(wt_data, mut_data)
    _write_results(results, outfile)


def _write_results(results, outfile):
    with open(outfile, 'w') as fh:
        fh.write('label, pvalue, tscore\n')
        for r in results:
            fh.write("{},{},{}\n".format(r[0], r[1], r[2]))

def _tstats(wt_data, mut_data):

    results = []

    for label, wt_values in wt_data.iteritems():
        mut_values = mut_data.get(label)
        if not mut_values:
            results.append((label, 'no mutant data', 'no mutant data'))
            continue

        tscore, pvalue = stats.ttest_ind(np.array(
            mut_values, dtype=np.float),
            np.array(wt_values, dtype=np.float))

        results.append((label, pvalue, tscore))

    # Sort results. Lowest pval first
    sorted_results = reversed(sorted(results, key=lambda pval: pval[2]))
    return sorted_results




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
    parser.add_argument('-o', dest='out', help='out file', required=True)
    args = parser.parse_args()

    organvolume_stats(args.wt, args.mut, args.out)