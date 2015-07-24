
"""
Convert the points section (need to manually extract that) from a vtk and spits out points that can be used
in elastix transformation. Need to manually add to top of output file:

index
numpoints


"""


in_ = '/home/neil/work/fake-abnormalities/vtk_files/right_thyroid.points'
out = '/home/neil/work/fake-abnormalities/vtk_files/right_thyroid.elxpoints'

numpoints = 0
def points(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

with open(in_, 'r') as fin, open(out, 'w') as fout:
    for line in fin.readlines():
        l = line.split()
        for p in points(l, 3):
            numpoints += 1
            fout.write(" ".join(p))
            fout.write('\n')
print 'points', numpoints
#elastix ... -fp fixedPointSet.txt -mp movingPointSet.txt