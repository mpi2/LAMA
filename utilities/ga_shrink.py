#!/usr/bin/env python

"""
Neil: Reimplementing the deap ga algorithm as I think I c might be able to make it go quicker
"""


import multiprocessing
import numpy as np
import random
import SimpleITK as sitk
import os
import copy
import math
try:
    import matplotlib.pyplot as plt
except ImportError:
    usemat = False
else:
    usemat = True


class GaShrink(object):
    def __init__(self, label, jac_value,out_dir,ngen, tourn_size=10, mutate_prob=0.0001, cross_over_chunk_size=14,
                 popsize=100, rand_mut=0.15):
        l = sitk.ReadImage(label)

        self.ideal = sitk.GetArrayFromImage(l).astype(np.float32)
        jac_float = float(jac_value)
        self.ideal[self.ideal == 1] = jac_float
        self.ideal[self.ideal == 0] = 1.0

        self.def_shape = list(self.ideal.shape) + [3]

        self.cross_over_chunk_size = cross_over_chunk_size

        self. mutate_prob = mutate_prob

        self.tourn_size = tourn_size

        self.popsize = popsize

        self.rand_mut = rand_mut
        self.run(out_dir, int(ngen))

    def evaluate(self, population):
        fits = []
        for ind in population:
            fits.append(self.calc_fitness(ind))
        return fits

    def calc_fitness(self, individual):
        jac_array = self.make_jac(individual)
        comp = np.sum(np.square(self.ideal - jac_array)) / individual.size
        return comp

    def make_jac(self, individual):

        im = sitk.GetImageFromArray(individual)
        jac = sitk.DisplacementFieldJacobianDeterminant(im)
        jac_array = sitk.GetArrayFromImage(jac)
        return jac_array

    def cross_over(self, ind1, ind2):
        """
        Given two indiduals (numpy arrays, switch differnt chunks over
        """

        r = np.random.randint(1, self.cross_over_chunk_size)
        numslices = ind1.shape[0]
        chunksize = int(math.floor(numslices / r))

        ind1[0:chunksize] = ind2[0:chunksize]
        ind2[chunksize:] = ind1[chunksize:]

        return ind1, ind2

    def mutate(self, ind):

        # Get a list of random vectors to mutate

        num_components = ind.size
        num_indices_to_mutate = int(num_components * self.mutate_prob)

        vector_indices = random.sample(xrange(0, num_components), num_indices_to_mutate)  # Mutate at these indices
        unraveled_vectors = ind.ravel()
        unraveled_vectors[vector_indices] += self.get_rand_mut_num()

        return ind


    def get_rand_mut_num(self):
        r = random.uniform(-self.rand_mut, self.rand_mut)
        return r


    def vector_to_def(self, vector):
        reshaped_ind = np.array(vector).reshape(self.def_shape)
        return reshaped_ind


    def get_initial_population(self, popsize, shape, val):
        pop = []
        for i in range(popsize):
            pop.append(np.random.uniform(-val, val, shape).astype(np.float32))
        return pop

    def generation_maker(self, pop, fits):

        #pick the top 10 to keep
        new_gen = []
        new_gen.extend([np.copy(x) for x in pop[0:5]])
        m = 5  # every 10

        p = 5
        while p < len(pop):

            if p %m != 0:

                fit1 = self.pick_fit(pop, fits)
                fit2 = self.pick_fit(pop, fits)
                prog1, prog2 = self.cross_over(np.copy(fit1), np.copy(fit2))
                new_gen.append(prog1)
                new_gen.append(prog2)
                p += 2
            else:
                if p+2 < len(pop):
                    progeny = np.copy(pop[p])
                    mut = self.mutate(progeny)
                    new_gen.append(mut)  # Make a copy of the indivual (breeding)
                    p += 1
        return new_gen


    def pick_fit(self, pop, fits):
        r_indexes = random.sample(range(len(fits)), self.tourn_size)
        r = min(r_indexes)
        return pop[r]

    def run(self, out_dir, ngen):

        #Prinrt the ideal jac
        ideal_out = os.path.join(out_dir, 'ideal_jac.nrrd')
        ideal_im = sitk.GetImageFromArray(self.ideal)
        sitk.WriteImage(ideal_im, ideal_out)
        ideal_shape = self.ideal.shape
        vector_field_shape = list(ideal_shape)
        vector_field_shape.append(3)  # 2 for 2d vectors

        # A list of np arrays. Contains individuals. These are numpy arrays each the size of the ideal * vector size
        population = self.get_initial_population(self.popsize, vector_field_shape, 0.1)

        temp_results = os.path.join(out_dir, 'intermediate_results')
        if not os.path.exists(temp_results):
            os.mkdir(temp_results)

        intermediate_num = 50  # Spit out intermediate results every x num generations

        with open(os.path.join(out_dir, 'chart_data.txt'), 'w')as dh:

            for gen in range(ngen):
                #print gen

                fits = self.evaluate(population)

                ordered_population = [x for (y, x) in sorted(zip(fits, population), key=lambda pair: pair[0])]
                fits = [y for (y, x) in sorted(zip(fits, population), key=lambda pair: pair[0])]
                population = self.generation_maker(ordered_population, fits)

                # population = ordered_population[0:popsize]

                dh.write('{}\n'.format(fits[0]))
                dh.flush()

                if gen != 0 and gen % intermediate_num == 0:
                    top10 = population[0:10]
                    self.write_results(temp_results, top10)
            top10 = population[0:10]
            self.write_results(out_dir, top10)

    def write_results(self, out, top10):
        for i, final in enumerate(top10):
            jac_result = self.make_jac(final)

            # write the jac
            outpath = os.path.join(out, "jac_{}.nrrd".format(i))
            sitk.WriteImage(sitk.GetImageFromArray(jac_result), outpath)

            # write the def
            outpath = os.path.join(out, "def_{}.nrrd".format(i))
            sitk.WriteImage(sitk.GetImageFromArray(self.vector_to_def(final)), outpath)





if __name__ == '__main__':
    import sys

    label = sys.argv[1]
    jac_value = sys.argv[2]
    out_dir = sys.argv[3]
    ngen = sys.argv[4]

    GaShrink(label, jac_value, out_dir, ngen)


