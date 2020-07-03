




import argparse
raise SystemExit('This CLI interafce needs updating')

parser = argparse.ArgumentParser()
parser.add_argument('-w', dest='wt', help='wt dir', required=True)
parser.add_argument('-m', dest='mut', help='mut dir', required=True)
parser.add_argument('-o', dest='output', help='output dir', required=True)
parser.add_argument('-s', dest='starts', help='start indices (x, y, z)', required=True, nargs=3, type=int)
parser.add_argument('-e', dest='ends', help='end indices (x, y, z)', required=True, nargs=3, type=int)

args = parser.parse_args()
wt = common.get_file_paths(args.wt)

mut = common.get_file_paths(args.mut)
normalise(wt, mut, args.elx_points, args.starts, args.ends)
