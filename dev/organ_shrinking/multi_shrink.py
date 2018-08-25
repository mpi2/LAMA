from argparse import ArgumentParser
from . import hill_climber_c
from multiprocessing import Pool
from os import listdir, mkdir, makedirs
from os.path import join, isdir
from glob import glob
import numpy as np


def run(volumes, organs, jacobian_range, output_folder, padding=5):

    # Loop through volumes
    jobs = []
    for folder in listdir(volumes):

        label_map = glob(join(volumes, folder, 'seg_*.nrrd'))[0]

        for label_num in organs:
            for jac_value in np.arange(jacobian_range[0], jacobian_range[1], 0.1, dtype=float):

                out_dir = join(output_folder, folder, 'jac{}'.format(int(jac_value*100)), 'label{}'.format(label_num))
                if not isdir(out_dir):
                    makedirs(out_dir)

                jobs.append({'label_map': label_map, 'label_num': label_num, 'jac_value': jac_value,
                             'padding': padding, 'out_dir': out_dir, 'ngen': 1000000, 'folder': folder})

    pool = Pool(4)
    list(zip(*pool.map(climb, jobs)))


def climb(job):

    print("Volume: {folder}\nLabel: {label_num}\nJacobian: {jac_value}".format(**job))
    hill_climber_c.HcShrink(job['label_map'], job['label_num'], job['jac_value'], job['padding'], job['out_dir'], job['ngen'])

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('-v', '--volumes', dest='volumes', required=True, help="Folder of embryos/segmentations")
    parser.add_argument('-l', '--labels', dest='labels', required=True, help="Organ labels to shrink", nargs='+')
    parser.add_argument('-j', '--jacobian', dest='jacobian', type=float, required=True,
                        help="Upper and lower Jacobian values to test", nargs=2)
    parser.add_argument('-o', '--output', dest='output', required=True, help="Output folder")

    args = parser.parse_args()
    run(args.volumes, args.labels, args.jacobian, args.output)



