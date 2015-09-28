#!/usr/bin/env python

from ga_shrink import GaShrink
import os
import subprocess



if __name__ == '__main__':
    label = '/home/neil/work/defined_abnormalites/270915_single_thymus/shrink_0.7_270915/thymus_label_roi.nrrd'
    jac_value = 0.7
    out_dir = '/home/neil/work/defined_abnormalites/270915_single_thymus/shrink_0.7_270915/optimize'
    out_file = '/home/neil/work/defined_abnormalites/270915_single_thymus/shrink_0.7_270915/optimize/optimize_results.txt'
    ngen = 400


    def test_parameter(name, options):

        rand_mut_res = []
        for m in options:

            GaShrink(label, jac_value, out_dir, ngen, **{name: m})

            data = os.path.join(out_dir, 'chart_data.txt')
            firstline = subprocess.check_output(['head', '-1', data])
            lastline = subprocess.check_output(['tail', '-1', data])

            f = float(firstline.strip())
            l = float(lastline.strip())

            rand_mut_res.append(f - l)
        print name, options, rand_mut_res
        return name, options, rand_mut_res

    with open(out_file, 'w')as fh:
        all_results = []

        res = test_parameter('cross_over_chunk_size', [11, 14, 18, 20])
        all_results.append(res)
        fh.write(str(all_results))
        fh.write('\n\n\n')
        fh.flush()
        print res



