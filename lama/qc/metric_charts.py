#!/usr/bin/env python

"""
Given a directory containing a multi-level registration, make plots of the cost function for each resolution.
Creates a html report concatentating all the images into one.
"""

ITER_PREFIX = 'IterationInfo'
CSS_FILE = 'style.css'

import os
from os.path import join, relpath
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from logzero import logger as logging


def make_charts(in_dir, out_dir):

    html = '<html><link rel="stylesheet" type="text/css" href="{}" /><body>\n'.format(CSS_FILE)
    stage = os.path.basename(in_dir)
    html += '<div class="title"> Elastix metric results for {}</div>\n'.format(stage)
    img_dir = join(out_dir, 'images')
    if not os.path.exists(img_dir):
        os.mkdir(img_dir)
    for subdir in os.listdir(in_dir):
        if not os.path.isdir(join(in_dir, subdir)):
            continue
        base = os.path.basename(subdir)
        html += '<div class="specimen clear">{}</div><br><br>'.format(base)

        files = [x for x in os.listdir(join(in_dir, subdir)) if x.startswith(ITER_PREFIX)]
        # Sort the iteration info files based on the resolution number preceding the .txt extension
        sorted_iter_files = sorted(files, key=lambda x: x.strip('.txt')[-1])
        for file_ in sorted_iter_files:
            if not file_.startswith(ITER_PREFIX):
                continue
            iter_path = join(in_dir, subdir, file_)
            iter_num = file_.split('.')[1]
            res_num = file_.split('.')[2][1]
            out_path = join(img_dir, base + '_' + iter_num + '_' + res_num + '.png')

            plot(iter_path, out_path)
            html += '<div class="chart"><img class="chart_img" src="{}"/></div>\n'.format(relpath(out_path, out_dir))
    html += '</body></html>'

    outfile = join(out_dir, 'iteration_report.html')
    with open(outfile, 'w') as html_fh:
        html_fh.write(html)
        html_fh.write(get_css())


def get_css():
    css = '<style>'
    css += 'body{font-family: Arial}'
    css += '.title{width: 100%; padding: 20px; background-color: lightblue; margin-bottom: 20px}'
    css += '.chart_img{width: 400px}\n'
    css += '.chart{float: left}'
    css += '.clear{clear: both}'
    css += '</style>'
    return css


def plot(iteration_file_, out_path):
    data = []

    with open(iteration_file_) as fh:
        # Get rid of header
        try:
            next(fh)
        except StopIteration as e:
            logging.warn(f'Problem reading iteration info from {iteration_file_} ')
            return
        for line in fh:
            metric = line.split()[1].strip()
            data.append(float(metric))
    plt.plot(data)
    plt.savefig(out_path)
    plt.close()


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser("Create metric charts for a stage of the registration")
    parser.add_argument('-i', '--in', dest='in_dir', help='Directory containing a series of registrations', required=True)
    parser.add_argument('-o', '--out', dest='out_file', help='File to export charts to', required=True)
    args = parser.parse_args()
    make_charts(args.in_dir, args.out_file)