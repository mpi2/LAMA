import nose
from nose import with_setup
from run_lama.stats import automated_annotation
import numpy as np
import tempfile
from os.path import join, abspath, dirname
import SimpleITK as sitk
import numpy as np
from os.path import join, realpath, dirname, isdir
import pandas as pd


def make_labels():
    """
    Make a test label file and write to a temporary nrrd. 9 labels
    """
    np.random.seed(44)
    # Make a volume with 5 labels 1-5
    labels = np.zeros(shape=(50, 50, 50))
    stats = np.zeros(shape=(50, 50, 50))
    for i in range(5):
        labels[i:i + 10, i:i + 10, i:i + 10] = i + 1
        stats[i:i + 10, i:i + 10, i:i + 10] = np.random.random_sample(1000).reshape((10,10,10))





    return labels, stats


def test_with_csv_file():
    """

    """

    label_info_file = tempfile.NamedTemporaryFile()
    label_info = """label,label_name,term
1,l1,emap:1
2,l2,emap:2
3,l3,emap:3
4,l4,emap:4
5,l5,emap:5
"""


    label_map, stats = make_labels()

    # Write out the label info file
    with open(label_info_file.name, 'w') as lif:
        lif.write(label_info)
    label_names = pd.read_csv(label_info_file.name)

    outfile = join(dirname(abspath(__file__)), 'test_data', 'test_output', 'autoannotator_simple.csv')

    # Here I'm using the labelmap as the dummy t-stats file as well. It should produce top hits for the
    #  higher -labelled organs
    ann = automated_annotation.Annotator(label_map, label_names, stats, outfile)
    ann.annotate()
