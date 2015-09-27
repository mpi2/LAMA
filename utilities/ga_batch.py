#!/usr/bin/env python

from ga_shrink import GaShrink
import os
import subprocess



if __name__ == '__main__':
    label = '/home/neil/work/defined_abnormalites_temp/230915_thymus/rewrite_ga/restart_250915/thymus_label_roi.nrrd'
    jac_value = 0.7
    out_dir = '/home/neil/work/defined_abnormalites_temp/230915_thymus/rewrite_ga/restart_250915/out'
    out_file = '/home/neil/work/defined_abnormalites_temp/230915_thymus/rewrite_ga/restart_250915/out/optimize_results.txt'
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


        res = test_parameter('mutate_prob', [0.001, 0.005, 0.0005, 0.0001])
        all_results.append(res)

        res = test_parameter('cross_over_chunk_size', [ 7, 8, 10, 14])
        all_results.append(res)
        fh.write(str(all_results))
        fh.write('\n\n\n')
        fh.flush()
        print res



