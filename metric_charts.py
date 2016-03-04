#!/usr/bin/env python


ITER_PREFIX = 'IterationInfo'
CSS_FILE = 'style.css'

import os
from os.path import join, relpath
import matplotlib.pyplot as plt

def make_charts(in_dir, out_dir):

    html = '<html><link rel="stylesheet" type="text/css" href="{}" /><body>\n'.format(CSS_FILE)
    img_dir = join(out_dir, 'images')
    if not os.path.exists(img_dir):
        os.mkdir(img_dir)
    for subdir in os.listdir(in_dir):
        if not os.path.isdir(join(in_dir, subdir)):
            continue
        base = os.path.basename(subdir)
        html += '<div class="specimen clear">{}</div><br><br>'.format(base)
        for file_ in os.listdir(join(in_dir, subdir)):
            if not file_.startswith(ITER_PREFIX):
                continue
            iter_path = join(in_dir, subdir, file_)
            iter_num = file_.split('.')[1]
            res_num = file_.split('.')[2][1]
            out_path = join(img_dir, base + '_' + iter_num + '_' + res_num + '.png')

            plot(iter_path, out_path)
            html += '<div class="chart"><img class="chart_img" src="{}"/></div>\n'.format(relpath(out_path, out_dir))
    html += '</body></html>'

    with open(join(out_dir, 'iteration_report.html'), 'w') as html_fh:
        html_fh.write(html)

    print_css(out_dir)

def print_css(out_dir):
    css = '.chart_img{width: 400px}\n'
    css += '.chart{float: left}'
    css += '.clear{clear: both}'
    # css += '.specimen{float: left}'
    with open(join(out_dir, CSS_FILE), 'w') as fh:
        fh.write(css)


def plot(iteration_file_, out_path):
    data = []

    with open(iteration_file_) as fh:
        # Get rid of header
        next(fh)
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