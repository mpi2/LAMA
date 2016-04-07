import struct
import numpy as np

def numpy_to_dat(mat, outfile):

    # create a binary file
    binfile = file(outfile, 'wb')
    # and write out two integers with the row and column dimension

    header = struct.pack('2I', mat.shape[0], mat.shape[1])
    binfile.write(header)
    # then loop over columns and write each
    for i in range(mat.shape[1]):
        data = struct.pack('%id' % mat.shape[0], *mat[:, i])
        binfile.write(data)

    binfile.close()

test_file = '/home/neil/work/circular_stats/test_cir.dat'

wt = []
mut = []
for i in range(12):
    wt.append(np.random.random_sample(10))
for i in range(4):
    mut.append(np.random.random_sample(10))

mat = np.vstack((wt, mut))

numpy_to_dat(mat, test_file)