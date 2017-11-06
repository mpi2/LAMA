#!/usr/bin/env python

"""
automated_annotation.py

Takes a label map, label names (and ontology) file and a 3D heatmap file (FDR filterd t-statistics for example) 
produces a csv of regions that are significantly affected accoring the statistics file.

Currently uses the following formula to create a score for each label (for both negative and positve hits to the same organ)

    neg_ratio = float(len(neg_t)) / float(label_vol)
    neg_score = median_neg_t * neg_ratio
"""
from argparse import ArgumentParser
import numpy as np
import math
import pandas as pd
import SimpleITK as sitk


class Annotator(object):

    def __init__(self, label_map_path, label_info_path, stats_path, outpath, type='jacobians'):
        """

        Parameters
        ----------
        label_map: path to label map
        label_info: path to label names file
            columns = label_number, label_name, emapa term (for example)
        stats: numpy ndarry
            FDR-thresholded t-statistics
        mask:
        """
        self.stats = path_to_array(stats_path)
        self.label_info = common.load_label_map_names(label_info_path)
        self.label_map = common.img_path_to_array(label_map_path)
        self.type = type

        self.no_labels = len(self.labelmap)

    def annotate(self):

        annotations = []

        for index, row in self.label_info.iterrows():
            label_num = row['label_num']
            description = row['description']
            term = row['term']

            label_mask = self.labelmap == label_num
            stats_organ = self.stats[label_mask]
            # Leave this out for now. James did this for organs that don't have any sidedness information in their name
            # side = 'right' if organ['right_label'] == str(label) else 'left'
            # side = '' if organ['right_label'] == organ['left_label'] else side

            # Calculate mean positive/negative t-statistics
            neg_t = stats_organ[stats_organ < 0]
            pos_t = stats_organ[stats_organ > 0]

            median_neg_t = np.nanmax((0, np.abs(np.median(neg_t))))
            median_pos_t = np.nanmax((0, np.abs(np.median(pos_t))))

            label_vol = np.sum(label_mask)

            if median_neg_t == 0:
               neg_score = 0
            else:
                neg_ratio = float(len(neg_t)) / float(label_vol)
                neg_score = median_neg_t * neg_ratio

            if median_pos_t == 0:
                pos_ratio = 0
                pos_score = 0
            else:
                pos_ratio = float(len(pos_t)) / float(label_vol)
                pos_score = median_pos_t * pos_ratio

            score = max(pos_score, neg_score)
            annotations.append({'label': label_num, 'name': description, 'term': term
                                'ratio': pos_ratio,
                                'median_pos_t': median_pos_t, 'median_neg_t': median_neg_t,
                                'score': score})

        # Create sorted dataframe
        df = pd.DataFrame(data=annotations).set_index('label').sort_values(['score'], ascending=False)
        df.to_csv(self.out_path)


def path_to_array(path):
    return sitk.GetArrayFromImage(sitk.ReadImage(path))

if __name__ == "__main__":
    import sys
    from os.path import join, dirname
    sys.path.insert(0, join(dirname(__file__), '..'))
    import common

    parser = ArgumentParser()
    parser.add_argument('-l', '--labelmap', dest='labelmap', help="Labelmap volume", required=True)
    parser.add_argument('-n', '--labelnames', dest='labelnames', help="CSV label names", required=True)
    parser.add_argument('-s', '--statistics', dest='stats', help="T-statistic volume", required=True)
    parser.add_argument('-o', '--outpath', dest='outpath', help="Path to save CSV to", default=None, required=True)
    args = parser.parse_args()

    ann = Annotator(args.labelmap,
                    args.label_names,
                    args.stats,
                    args.outpath)

    df = ann.annotate()
    print(df)
