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
import pandas as pd
import SimpleITK as sitk


class Annotator(object):

    def __init__(self, label_map: np.ndarray, label_info, stats, outpath, type='jacobians'):
        """

        Parameters
        ----------
        label_map: numpy.ndarray labelmap
        label_info: pandas dataframe
            example
                label,label_name,term
                1,l1,emap:1
                2,l2,emap:2
                3,l3,emap:4
                4,l4,emap:5
                5,l5,emap:6
        stats: numpy ndarry
            FDR-thresholded t-statistics volume
        outpath: str
            path to outfile
        mask:
        """
        self.stats = stats
        self.label_info = label_info
        self.labelmap = label_map
        self.type = type
        self.outpath = outpath
        if label_map.shape != stats.shape:
            raise ValueError("Annotator: label map shape {} does not match heatmap shape {}".format(
                label_map.shape, stats.shape
            ))

    def annotate(self):

        annotations = []  # each annotation saved here as a dict and made into a pandas dataframe at the end

        for i, v in self.label_info.iterrows():
            label_num = int(v['label'])
            description = v['label_name']
            term = v['term']

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
            annotations.append({'label': label_num, 'name': description, 'term': term,
                                'ratio': pos_ratio,
                                'median_pos_t': median_pos_t, 'median_neg_t': median_neg_t,
                                'score': score})

        # Create sorted dataframe
        df = pd.DataFrame(data=annotations).set_index('label').sort_values(['score'], ascending=False)
        df.to_csv(self.outpath)
        return df


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

    ann = Annotator(path_to_array(args.labelmap),
                    common.load_label_map_names(args.labelnames, include_terms=True),
                    path_to_array(args.stats),
                    args.outpath)

    df = ann.annotate()
    print(df)
