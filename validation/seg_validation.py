from __future__ import division
from argparse import ArgumentParser
from os.path import basename, join, isfile
import common
import SimpleITK as sitk
import csv
from scipy.spatial.distance import jaccard
from scipy.spatial.distance import dice
from scipy.spatial.distance import rogerstanimoto
import pandas as pd
import seaborn as sns


class LabelValidation(object):

    def __init__(self, inv_dir, ground_truth, out_dir, organ_names):

        self.inverted_volumes = common.GetFilePaths(inv_dir)
        self.ground_truth_volumes = common.GetFilePaths(ground_truth)
        self.out_dir = out_dir

        with open(organ_names, 'rb') as org:
            self.organs = [row for row in csv.DictReader(org)]

        self.score_data = join(self.out_dir, 'scores.csv')
        if isfile(self.score_data):
            with open(self.score_data, 'rb') as f:
                scores = pd.DataFrame.from_csv(f)
        else:
            scores = self.run()

        # Plot results
        self.plot_results(scores)

    def run(self):

        scores = pd.DataFrame()
        for true_path in self.ground_truth_volumes:

            inv_path = self.find_matching(true_path)
            if not inv_path:
                continue

            vol_name = basename(inv_path)
            print "Loading volumes: '{}'".format(vol_name)

            # Run label statistics
            true_im, inv_im = sitk.ReadImage(true_path), sitk.ReadImage(inv_path)
            true_stats, inv_stats = sitk.LabelStatisticsImageFilter(), sitk.LabelStatisticsImageFilter()
            true_stats.Execute(true_im, true_im)
            inv_stats.Execute(inv_im, inv_im)

            # Convert to np arrays
            true_arr = sitk.GetArrayFromImage(true_im).flatten()
            inv_arr = sitk.GetArrayFromImage(inv_im).flatten()

            for row in self.organs:

                # Get labels to compare
                organ_name = row['name']
                true_label, inv_label = int(row['true_label']), int(row['inv_label'])

                # Scores
                ratio_score = true_stats.GetCount(true_label) / inv_stats.GetCount(inv_label)
                jaccard_score = 1 - jaccard(true_arr == true_label, inv_arr == inv_label)
                dice_score = 1 - dice(true_arr == true_label, inv_arr == inv_label)
                tanimoto_score = 1 - rogerstanimoto(true_arr == true_label, inv_arr == inv_label)

                scores = scores.append({'organ_name': organ_name, 'ratio': ratio_score, 'jaccard': jaccard_score,
                                        'dice': dice_score, 'tanimoto': tanimoto_score}, ignore_index=True)

        # Dump and return results
        with open(self.score_data, 'wb') as f:
            scores.to_csv(f)
        return scores

    def plot_results(self, scores):

        for metric in ['dice', 'jaccard', 'ratio', 'tanimoto']:
            sns.swarmplot(x='organ_name', y=metric, data=scores)
            sns.plt.title(metric.title())
            sns.plt.savefig(join(self.out_dir, "{}.png".format(metric)), dpi=600)
            sns.plt.clf()

    def find_matching(self, g_vol):

        try:
            return [vol for vol in self.inverted_volumes if basename(g_vol) in vol][0]
        except IndexError:
            print "No inverted volume found for '{}'".format(basename(g_vol))
            return None


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument('-i', '--inverted', dest='inv_dir', required=True, help="Directory of inverted labelmaps")
    parser.add_argument('-g', '--ground-truth', dest='ground_truth', required=True, help="Directory of manual ground truth")
    parser.add_argument('-n', '--names', dest='organ_names', required=True, help="CSV file of organ names")
    parser.add_argument('-o', '--output', dest='out_dir', required=True, help="Directory to output validation results")
    args = parser.parse_args()

    val = LabelValidation(args.inv_dir, args.ground_truth, args.out_dir, args.organ_names)