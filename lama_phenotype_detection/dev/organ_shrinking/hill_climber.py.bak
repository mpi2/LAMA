#!/usr/bin/env python

"""
Henrik: Reimplementing this to be a hill climber
"""


import numpy as np
import SimpleITK as sitk
import os
import copy


class HcShrink(object):
    def __init__(self, label, jac_value, out_dir,ngen):
        l = sitk.ReadImage(label)

        self.ideal = sitk.GetArrayFromImage(l).astype(np.float32)
        jac_float = float(jac_value)
        self.ideal[self.ideal == 1] = jac_float
        self.ideal[self.ideal == 0] = 1.0

        # Henrik: not sure what this is for
        self.def_shape = list(self.ideal.shape) + [3]

        self.run(out_dir, int(ngen))

    @staticmethod
    # get the determinents per voxel
    def make_jac(individual):

        im = sitk.GetImageFromArray(individual)
        jac = sitk.DisplacementFieldJacobianDeterminant(im)
        jac_array = sitk.GetArrayFromImage(jac)
        return jac_array

    # calculate the fitness
    def calc_fitness(self, individual):
        """
        Add weight to regions where ideal_jac != 1.0
        """
        jac_array = self.make_jac(individual)
        surronding_comp = np.sum(np.square(self.ideal[self.ideal == 1.0] - jac_array[self.ideal == 1.0])) / individual.size
        organ_comp = np.sum(np.square(self.ideal[self.ideal != 1.0] - jac_array[self.ideal != 1.0])) / individual.size
        return surronding_comp + (organ_comp *3)

    # evalualte a whole population
    def evaluate(self, ind):
        return self.calc_fitness(ind)

    def mutate(self, ind, fit):

        # pick a random position
        random_position = np.random.randint(0, ind.size)
        mutate_val = np.random.uniform(0, fit *  30)
        rand = np.random.uniform(-mutate_val, mutate_val)
        ind.ravel()[random_position] += rand
        return ind

    def vector_to_def(self, vector):
        reshaped_ind = np.array(vector).reshape(self.def_shape)
        return reshaped_ind

    def get_initial_individual(self, shape, val):

        ind = np.zeros(shape, dtype=np.float32)
        #ind = np.random.uniform(-val, val, shape).astype(np.float32)
        #ind[self.ideal == 1.0] = 0.0
        return ind

    def run(self, out_dir, ngen):

        #Print the ideal jac
        ideal_out = os.path.join(out_dir, 'ideal_jac.nrrd')
        ideal_im = sitk.GetImageFromArray(self.ideal)
        sitk.WriteImage(ideal_im, ideal_out)
        ideal_shape = self.ideal.shape

        vector_field_shape = list(ideal_shape)
        vector_field_shape.append(3)  # 2 for 2d vectors

        print vector_field_shape

        # A list of np arrays. Contains individuals. These are numpy arrays each the size of the ideal * vector size
        ind = self.get_initial_individual(vector_field_shape, 0.0)

        print ind.shape

        temp_results = os.path.join(out_dir, 'intermediate_results')
        if not os.path.exists(temp_results):
            os.mkdir(temp_results)

        intermediate_num = 10  # Spit out intermediate results every x num generations

        with open(os.path.join(out_dir, 'chart_data.txt'), 'w')as dh:

            for gen in range(ngen):

                # evaluate individual
                fit = self.evaluate(ind)

                #print fit

                temp = True
                # mutate the individual until it gets better
                while temp:
                    temp_ind = copy.deepcopy(ind)
                    temp_ind = self.mutate(temp_ind, fit)

                    if self.evaluate(temp_ind) < fit:
                        ind = temp_ind
                        temp = False
                dh.write('{}\n'.format(fit))
                dh.flush()

                if gen != 0 and gen % intermediate_num == 0:
                    self.write_results(temp_results, ind)
            self.write_results(out_dir, ind)

    def write_results(self, out, ind):

        jac_result = self.make_jac(ind)

        # write the jac
        outpath = os.path.join(out, "jac_{}.nrrd".format(1))
        sitk.WriteImage(sitk.GetImageFromArray(jac_result), outpath)

        # write the def
        outpath = os.path.join(out, "def_{}.nrrd".format(1))
        sitk.WriteImage(sitk.GetImageFromArray(self.vector_to_def(ind)), outpath)





if __name__ == '__main__':
    import sys

    label = sys.argv[1]
    jac_value = sys.argv[2]
    out_dir = sys.argv[3]
    ngen = sys.argv[4]

    HcShrink(label, jac_value, out_dir, ngen)

