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


def evaluate(population):
    fits = []
    for ind in population:
        fits.append(calc_fitness(ind))
    return fits


def calc_fitness(individual):
    jac_array = make_jac(individual)
    comp = np.sum(np.square(ideal - jac_array)) / individual.size
    return comp

def make_jac(individual):

    im = sitk.GetImageFromArray(individual)
    jac = sitk.DisplacementFieldJacobianDeterminant(im)
    jac_array = sitk.GetArrayFromImage(jac)
    return jac_array

def cross_over(ind1, ind2):
    """
    Given two indiduals (numpy arrays, switch differnt chunks over
    """
    r = np.random.randint(1, 4)
    numslices = ind1.shape[0]
    chunksize = int(math.floor(numslices / r))

    ind1[0:chunksize] = ind2[0:chunksize]
    ind2[chunksize:] = ind1[chunksize:]


def mutate(ind, indpb=0.01):

    # Get a list of random vectors to mutate
    shape = ind.shape
    num_vectors = (np.prod(ind.shape[:-1]))  # Number of vectors to mutate
    num_indices_to_mutate = int(num_vectors * indpb)

    vector_indices = random.sample(xrange(0, num_vectors), num_indices_to_mutate)  # Mutate at these indices

    # add weight to either changing magnitude or direction
    magnitude_weight = 0.9
    r = np.random.uniform(0, 1)
    reshapaed = ind.reshape(np.prod(ind.shape[:-1]), 3)
    if r < magnitude_weight:

        reshapaed[vector_indices] += get_rand_mut_num()
    else:

        for indx in [vector_indices]:
                reshapaed[indx] += get_rand_mut_num()

    return ind


def get_rand_mut_num():
    rand_mut = 0.2
    r = random.uniform(-rand_mut, rand_mut)
    return r


def vector_to_def(vector):
    reshaped_ind = np.array(vector).reshape(def_shape)
    return reshaped_ind


def get_initial_population(popsize, shape, val):
    pop = []
    for i in range(popsize):
        pop.append(np.random.uniform(-val, val, shape).astype(np.float32))
    return pop


def generation_maker(pop):
    mutated = []

    breeding = 2  # breed every 10
    for i, ind in enumerate(pop):
        if i % breeding == 0:
            if i+2 < len(pop):
                cross_over(pop[i], pop[i+1])
        else:
            progeny = np.copy(ind)
            mutated.append(mutate(progeny))  # Make a copy of the indivual (breeding)
    pop.extend(mutated)




def run(out_dir, ngen):


    popsize = 80
    #Prinrt the ideal jac
    ideal_out = os.path.join(out_dir, 'ideal_jac.nrrd')
    ideal_im = sitk.GetImageFromArray(ideal)
    sitk.WriteImage(ideal_im, ideal_out)
    ideal_shape = ideal.shape
    vector_field_shape = list(ideal_shape)
    vector_field_shape.append(3)  # 2 for 2d vectors

    # A list of np arrays. Contains individuals. These are numpy arrays each the size of the ideal * vector size
    population = get_initial_population(popsize, vector_field_shape, 0.2)

    temp_results = os.path.join(out_dir, 'intermediate_results')
    if not os.path.exists(temp_results):
        os.mkdir(temp_results)

    intermediate_num = 50  # Spit out intermediate results every x num generations

    with open(os.path.join(out_dir, 'chart_data.txt'), 'w')as dh:

        for gen in range(ngen):
            print gen
            generation_maker(population)
            fits = evaluate(population)

            ordered_population = [x for (y, x) in sorted(zip(fits, population), key=lambda pair: pair[0])]

            population = ordered_population[0:popsize]

            dh.write('{}\n'.format(fits[0]))
            dh.flush()

            if gen != 0 and gen % intermediate_num == 0:
                 top10 = population[0:10]
                 write_results(temp_results, top10)

        write_results(out_dir, top10)

def write_results(out, top10):
    for i, final in enumerate(top10):
        jac_result = make_jac(final)

        # write the jac
        outpath = os.path.join(out, "jac_{}.nrrd".format(i))
        sitk.WriteImage(sitk.GetImageFromArray(jac_result), outpath)

        # write the def
        outpath = os.path.join(out, "def_{}.nrrd".format(i))
        sitk.WriteImage(sitk.GetImageFromArray(vector_to_def(final)), outpath)

if __name__ == '__main__':
    import sys

    label = sys.argv[1]
    jac_value = sys.argv[2]
    out_dir = sys.argv[3]
    ngen = sys.argv[4]

    l = sitk.ReadImage(label)
    ideal = sitk.GetArrayFromImage(l).astype(np.float32)
    jac_float = float(jac_value)
    ideal[ideal == 1] = jac_float
    ideal[ideal == 0] = 1.0
    def_shape = list(ideal.shape) + [3]

    run(out_dir, int(ngen))

