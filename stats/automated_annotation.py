from argparse import ArgumentParser
import numpy as np
import math
import pandas as pd
import SimpleITK as sitk


class Annotator(object):

    def __init__(self, label_map, label_names, stats, outpath, mask=None):
        """

        Parameters
        ----------
        label_map: numpy ndarray
        label_names: dict
            {o: 'organ name', 1: 'organ name' ....}
        stats: numpy ndaarry
            FDR-thresholded t-statistics
        mask:
        """

        self.out_path = outpath
        self.labelmap = label_map
        self.stats = stats
        self.mask = mask
        self.label_names = label_names

        self.no_labels = len(self.labelmap)

    def annotate(self):

        annotations = []

        for organ_num, organ_name in self.label_names.iteritems():

            label_mask = self.labelmap == organ_num
            stats_organ = self.stats[label_mask]
            # Leave this out for now. James did this for organs that don't have any sidedness information in their name
            # side = 'right' if organ['right_label'] == str(label) else 'left'
            # side = '' if organ['right_label'] == organ['left_label'] else side

            # Calculate mean positive/negative t-statistics
            neg_t = stats_organ[stats_organ < 0]
            pos_t = stats_organ[stats_organ > 0]

            median_neg_t = np.abs(np.median(neg_t))
            median_pos_t = np.abs(np.median(pos_t))

            median_neg_t = 0 if math.isnan(median_neg_t) else median_neg_t
            median_pos_t = 0 if math.isnan(median_pos_t) else median_pos_t

            # Calculate volume ratios
            label_vol = np.sum(label_mask)
            neg_ratio = float(len(neg_t)) / float(label_vol)
            pos_ratio = float(len(pos_t)) / float(label_vol)

            # Calculate "score"
            neg_score = median_neg_t * neg_ratio
            pos_score = median_pos_t * pos_ratio

            annotations.append({'label': organ_num, 'name': organ_name,
                                'pos_ratio': pos_ratio, 'neg_ratio': neg_ratio,
                                'median_pos_t': median_pos_t, 'median_neg_t': median_neg_t,
                                'pos_score': pos_score, 'neg_score': neg_score})

        # Create sorted dataframe
        df = pd.DataFrame(data=annotations).set_index('label').sort_values(['pos_score', 'neg_score'], ascending=False)
        df.to_csv(self.out_path)
        return df


def path_to_array(path):
    return sitk.GetArrayFromImage(sitk.ReadImage(path))

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument('-l', '--labelmap', dest='labelmap', help="Labelmap .mnc", required=True)
    parser.add_argument('-n', '--labelnames', dest='labelnames', help="CSV label names", required=True)
    parser.add_argument('-s', '--statistics', dest='stats', help="T-statistic volume", required=True)
    parser.add_argument('-m', '--mask', dest='mask', help="Statistisc mask", default=None, required=False)
    parser.add_argument('-o', '--outpath', dest='outpath', help="Path to save CSV to", default=None, required=True)
    args = parser.parse_args()

    lns = {}
    with open(args.labelnames, 'rb') as fh:
        for i, line in enumerate(fh):
            lns[i +1] = line.strip()
    if args.mask:
        mask = path_to_array(args.mask)
    else:
        mask = None

    ann = Annotator(path_to_array(args.labelmap),
                    lns, path_to_array(args.stats),
                    outpath=args.outpath,
                    mask=mask)
    df = ann.annotate()
    print(df)
