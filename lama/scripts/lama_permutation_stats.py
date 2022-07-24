#! /usr/bin/env python3

"""
Usage
-----
This script should be run from the root directory (one up from this script)


todo: This bit
$ cd lama_phenotype_detection
$ scripts/lama_stats.py -c <path to stats config> -w <Path to wild type diretory>
"""
from pathlib import Path
import sys

import yaml

from lama.stats.permutation_stats import run_permutation_stats


allowed_cfg_keys = [
    'wildtype_dir',
    'mutant_dir',
    'output_dir',
    'treatment_dir',
    'interaction_dir',
    'n_permutations',
    'label_metadata',
    'label_map',
    'norm_to_whole_embryo_vol',
    'qc_file',
    'voxel_size',
    'two_way'
]


def main():
    import argparse

    parser = argparse.ArgumentParser("Permutation-based stats")
    parser.add_argument('-c', '--config', dest='cfg_path', help='wildtype registration directory', required=True,
                        type=str)

    args = parser.parse_args()
    run(args.cfg_path)


def run(cfg_path):

    def p(path):
        if path is None:
                return

        cfg_dir = Path(cfg_path).parent

        resolved = (cfg_dir / path).resolve()

        if not resolved.exists():
            raise FileNotFoundError(f'Cannot find: {resolved}')
        return resolved

    with open(cfg_path, 'r') as fh:
        cfg = yaml.load(fh)

    for key in cfg.keys():
        if key not in allowed_cfg_keys:
            raise ValueError(f'Config key "{key}" is not in the allowed keys: {", ".join(allowed_cfg_keys)} ')

    # required parameters
    try:
        wt_dir = p(cfg['wildtype_dir'])
        mut_dir = p(cfg['mutant_dir'])
    except KeyError:
        raise KeyError("'wildtype_dir', 'mutant_dir' are required parameters")

    out_dir = p(cfg.get('output_dir', Path(cfg_path).parent))

    # Optional parameters

    n_perm = int(cfg.get('n_permutations', 1000))
    label_meta = p(cfg.get('label_metadata'))
    label_map = p(cfg.get('label_map'))
    wev_norm = bool(cfg.get('norm_to_whole_embryo_vol', True))
    qc_file = p(cfg.get('qc_file'))
    voxel_size = float(cfg.get('voxel_size', 1.0))

    treat_dir = p(cfg['treatment_dir'])
    inter_dir = p(cfg['interaction_dir'])
    two_way = bool(cfg.get('two_way', False))

    run_permutation_stats.run(wt_dir=wt_dir,
                              mut_dir=mut_dir,
                              out_dir=out_dir,
                              num_perms=n_perm,
                              label_info=label_meta,
                              label_map_path=label_map,
                              normalise_to_whole_embryo=wev_norm, qc_file=qc_file,
                              voxel_size=voxel_size,
                              two_way=two_way,
                              treat_dir=treat_dir,
                              inter_dir=inter_dir
    )


if __name__ == '__main__':
    main()

