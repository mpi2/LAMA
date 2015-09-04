import sys

from mrch_regpipeline import invert_isosurfaces as iso

invert_config = sys.argv[1]
mesh_dir = sys.argv[2]
iso_out_dir = sys.argv[3]

iso(invert_config, mesh_dir, iso_out_dir)
