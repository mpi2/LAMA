#!/usr/bin/env python

"""
Adapted from here: https://github.com/Mouse-Imaging-Centre/minc-stuffs/blob/master/python/TFCE
"""


from optparse import OptionParser
from scipy import ndimage, weave
from numpy import *
import operator
import SimpleITK as sitk
import numpy as np


def tfce(invol, outvol, dh=0.1, E=0.5, H=2.0):
    # infile = "1.mnc"
    # out = "1-tfce.mnc"

    # invol = volumeFromFile(infile)
    # outvol = volumeFromInstance(invol, outfile, volumeType="ushort")

    # dh = 0.1
    # E = 0.5
    # H = 2.0
    # outvol.data  # make sure data has been loaded

    # step through data with dh increments
    for h in arange(0, invol.max(), dh):
        # threshold the data with current height
        thresh = array(invol > h, "uint8")

        # connected components labelling
        l = ndimage.label(thresh)
        print("L:", l[1])
        # compute the size of each label
        sizes = array(ndimage.sum(thresh, l[0], range(l[1] + 1)))
        # modulate label size by voxel volume
        #sizes = sizes * reduce(operator.mul, invol.separations)
        print("sizes", sizes.shape)

        # compute TFCE
        if l[1] > 0:
            print
            "inside", h, l[1]
            code = """
            for (int x=0; x < nx; ++x) {
               for (int y=0; y < ny; ++y) {
                  for (int z=0; z < nz; ++z) {
                     if (labeled(x,y,z) > 0) {
                        int e = labeled(x,y,z);
                        //std::cout << data(x,y,z) << " " << h << " " << sizes(e) << std::endl;
                        data(x,y,z) = data(x,y,z) + (pow(h, H) * pow(sizes(e), E) * dh);
                        //data(x,y,z) = n;
            }}}}
            """
            data = outvol
            nx, ny, nz = data.shape
            labeled = l[0]

            weave.inline(code, ['data', 'nx', 'ny', 'nz', 'labeled', 'h',
                                'H', 'E', 'sizes', 'dh'],
                         type_converters=weave.converters.blitz,
                         compiler='gcc')

        print(h)
        # outvol.writeFile()
        # outvol.closeVolume()


if __name__ == "__main__":

    usage = "usage: %prog [options] input.mnc output.mnc"
    description = """
Applies the texture free cluster enhancement (TFCE) to a statistics image. See:
Smith and Nichols. Threshold-free cluster enhancement: addressing problems of smoothing, threshold dependence and localisation in cluster inference. NeuroImage (2009) vol. 44 (1) pp. 83-98
    """
    parser = OptionParser(usage=usage, description=description)

    parser.add_option("-d", "--dh", dest="dh",
                      help="Increments over which to compute TFCE [default: %default]",
                      type="float", default=0.1)
    parser.add_option("-E", dest="E",
                      help="Power by which to raise the extent [default: %default]",
                      type="float", default=0.5)
    parser.add_option("-H", dest="H",
                      help="Power by which to raise the height [default: %default]",
                      type="float", default=2.0)
    parser.add_option("--pos-and-neg", dest="pos_and_neg",
                      help="Use both positive and negative data in input [default]",
                      action="store_const", const="both", default="both")
    parser.add_option("--pos-only", dest="pos_and_neg",
                      help="Use only positive data in input",
                      action="store_const", const="pos")
    parser.add_option("--neg-only", dest="pos_and_neg",
                      help="Use only negative data in input",
                      action="store_const", const="neg")

    (options, args) = parser.parse_args()



    invol_path = '/home/neil/sig/LAMA_results/E14.5/compare_cbx2_male_female_variances/male_female_and_cbx2/all_inputs/full_run_14um/stats_male_female/jacobians_male_female_192_to_12/LinearModelR/__jacobians_male_female_192_to_12_genotype_LinearModelR_Tstats_genotype_stats_.nrrd'
    invol = sitk.GetArrayFromImage(sitk.ReadImage(invol_path))
    outvol_path = '/home/neil/sig/LAMA_results/E14.5/compare_cbx2_male_female_variances/male_female_and_cbx2/all_inputs/full_run_14um/stats_male_female/jacobians_male_female_192_to_12/LinearModelR/__jacobians_male_female_192_to_12_genotype_LinearModelR_Tstats_genotype_stats_TFCE_.nrrd'
    outvol_pos = np.zeros_like(invol, dtype=np.short)
    outvol_neg = np.zeros_like(invol, dtype=np.short)



    tfce(invol, outvol_pos, options.dh, options.E, options.H)
    print('pos vol', outvol_pos.min(), outvol_pos.max())
    invol *= -1

    tfce(invol, outvol_neg, options.dh, options.E, options.H)
    print('neg ', outvol_neg.min(), outvol_neg.max())

    outvol_neg *= -1
    print('neg ', outvol_neg.min(), outvol_neg.max())

    outvol_pos += outvol_neg
    outimg = sitk.GetImageFromArray(outvol_pos)
    sitk.WriteImage(outimg, outvol_path)