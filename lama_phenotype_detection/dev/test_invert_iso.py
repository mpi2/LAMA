#!/usr/bin/env python

import mrch_regpipeline as reg

def run(config, mesh_dir, out_dir):
    reg.invert_isosurfaces(config, mesh_dir, out_dir)


if __name__ == '__main__':
    import sys

    config = sys.argv[1]
    mesh_dir = sys.argv[2]
    out_dir = sys.argv[3]

    run(config, mesh_dir, out_dir)


