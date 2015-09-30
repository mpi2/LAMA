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
manager = multiprocessing.Manager()
metric_list = manager.list()

def process_ga(label, jac_value, out_dir, ngen=500, numrounds=100, numthreads=4, tourn_size=10, mutate_prob=0.0001,
               cross_over_chunk_size=14, popsize=100, rand_mut=0.15):
    """
    """

    num_consumers = numthreads
    results_q = multiprocessing.Queue(1000000000)

    chart_data_path = os.path.join(out_dir, 'chart_data.txt')
    chart_file = open(chart_data_path, 'w')

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    def_results_dir = os.path.join(out_dir, 'defs')
    if not os.path.exists(def_results_dir):
        os.mkdir(def_results_dir)
    jac_results_dir = os.path.join(out_dir, 'jacs')
    if not os.path.exists(jac_results_dir):
        os.mkdir(jac_results_dir)

    print 'Creating %d consumers' % num_consumers
    try:
        consumers = list()
        for a in range(num_consumers):
            p = GaShrink(label, jac_value, out_dir, results_q, str(a), def_results_dir, jac_results_dir, ngen=ngen, tourn_size=10, mutate_prob=0.0001, cross_over_chunk_size=14,
                     popsize=100, rand_mut=0.15, initial=True)
            consumers.append(p)
            p.start()

        for w in consumers:
            w.join()

        for n in range(numrounds):
            consumers = list()
            for a in range(num_consumers):
                p = GaShrink(label, jac_value, out_dir, results_q, str(a), def_results_dir, jac_results_dir, ngen=10, tourn_size=10, mutate_prob=0.0001, cross_over_chunk_size=14,
                         popsize=100, rand_mut=0.15, initial=False)
                consumers.append(p)
                p.start()

            for w in consumers:
                w.join()

            chart_file.write(str(min(metric_list)) + '\n')
            chart_file.flush()

    except KeyboardInterrupt:
        print "Caught KeyboardInterrupt, terminating workers"
        for w in consumers:
            w.terminate()
            w.join()
    chart_file.close()


class GaShrink(multiprocessing.Process):
    def __init__(self, label, jac_value, out_dir, results_q, id, def_dir, jac_dir, ngen=1000, tourn_size=10, mutate_prob=0.0001, cross_over_chunk_size=14,
                 popsize=100, rand_mut=0.15, initial=True):
        super(GaShrink, self).__init__()

        l = sitk.ReadImage(label)

        self.initial = initial

        self.id = id
        self.def_dir = def_dir
        self.jac_dir = jac_dir

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

        self.results_q = results_q

        self.out_dir = out_dir

        self.ngen = int(ngen)

    def get_result(self):
        return self.top10

    def evaluate(self, population):
        fits = []
        for ind in population:
            fits.append(self.calc_fitness(ind))
        return fits

    def calc_fitness(self, individual):
        """
        Add weight to regions where ideal_jac != 1.0
        """
        jac_array = self.make_jac(individual)
        surronding_comp = np.sum(np.square(self.ideal[self.ideal == 1.0] - jac_array[self.ideal == 1.0])) / individual.size
        organ_comp = np.sum(np.square(self.ideal[self.ideal != 1.0] - jac_array[self.ideal != 1.0])) / individual.size
        return surronding_comp + (organ_comp / 3)

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
        while p < self.popsize:

            if p %m != 0 or len(pop) < self.popsize:

                fit1 = self.pick_fit(pop, fits)
                fit2 = self.pick_fit(pop, fits)
                prog1, prog2 = self.cross_over(np.copy(fit1), np.copy(fit2))
                new_gen.append(prog1)
                new_gen.append(prog2)
                p += 2
            else:
                if p+2 < self.popsize:
                    progeny = np.copy(pop[p])
                    mut = self.mutate(progeny)
                    new_gen.append(mut)  # Make a copy of the indivual (breeding)
                    p += 1
        return new_gen

    def pick_fit(self, pop, fits):
        r_indexes = random.sample(range(len(fits)), self.tourn_size)
        r = min(r_indexes)
        return pop[r]

    def get_previous_round(self):
        pop = []
        previous_result_paths = os.listdir(self.def_dir)
        for p in previous_result_paths:
            fullpath = os.path.join(self.def_dir, p)
            array = sitk.GetArrayFromImage(sitk.ReadImage((fullpath)))
            pop.append(array)
        return pop

    def run(self):

        ideal_shape = self.ideal.shape
        vector_field_shape = list(ideal_shape)
        vector_field_shape.append(3)  # 2 for 2d vectors

        # A list of np arrays. Contains individuals. These are numpy arrays each the size of the ideal * vector size

        if self.initial:
            population = self.get_initial_population(self.popsize, vector_field_shape, 0.01)
        else:
            population = self.get_previous_round()

        for gen in range(self.ngen):
            fits = self.evaluate(population)

            ordered_population = [x for (y, x) in sorted(zip(fits, population), key=lambda pair: pair[0])]
            fits = [y for (y, x) in sorted(zip(fits, population), key=lambda pair: pair[0])]
            population = self.generation_maker(ordered_population, fits)


        for i, ind in enumerate(population[0:5]):
            # write the ja
            jac_result = self.make_jac(ind)
            jac_filename = os.path.join(self.jac_dir, "jac_{}{}.nrrd".format(i, self.id))
            sitk.WriteImage(sitk.GetImageFromArray(jac_result), jac_filename)
            def_filename = os.path.join(self.def_dir, self.id + '_id' + str(i) + 'def_.nrrd')
            sitk.WriteImage(sitk.GetImageFromArray(ind), def_filename)

        metric_list.append(fits[0])
        return


if __name__ == '__main__':
    import sys

    label = sys.argv[1]
    jac_value = sys.argv[2]
    out_dir = sys.argv[3]
    ngen = sys.argv[4]
    thread_rounds = sys.argv[5]
    num_threads = sys.argv[6]


    process_ga(label, jac_value, out_dir, int(ngen), int(thread_rounds), int(num_threads))


