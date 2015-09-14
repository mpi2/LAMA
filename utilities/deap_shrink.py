#!/usr/bin/env python


import multiprocessing
from deap import creator, base, tools, algorithms
import numpy as np
import random
import SimpleITK as sitk
import os
try:
    import matplotlib.pyplot as plt
except ImportError:
    usemat = False
else:
    usemat = True




def calc_fitness(individual):
    jac_array = make_jac(individual)
    comp = (1 - np.sum(np.square(ideal - jac_array))) / len(individual)
    return [comp]

def make_jac(individual):
    reshaped_ind = vector_to_def(individual)
    im = sitk.GetImageFromArray(reshaped_ind)
    jac = sitk.DisplacementFieldJacobianDeterminant(im)
    jac_array = sitk.GetArrayFromImage(jac)
    return jac_array

def cross_over(ind1, ind2):

    step = 3
    for i in range(0, len(ind1) - step, step):
        trip1 = ind1[i: i+ step]
        trip2 = ind2[i: i+step]
        rnum = random.random()
        if rnum < 0.1:
            ind1[i: i+step] = trip2
            ind2[i: i+step] = trip1
    return ind1, ind2


def mutate(ind, indpb=0.05):

    step = 3

    for i in range(0, len(ind) - step, step):
        rnum = random.random()
        if rnum < indpb:
            ind[i] = ind[i] + get_rand_mut_num()
            ind[i + 1] = ind[i + 1] + get_rand_mut_num()
            ind[i + 2] = ind[i + 2] + get_rand_mut_num()
    return ind,

    # size = len(individual)
    # for i in xrange(size):
    #     if random.random() < indpb:
    #         swap_indx = random.randint(0, size - 2)
    #         if swap_indx >= i:
    #             swap_indx += 1
    #         individual[i], individual[swap_indx] = \
    #             individual[swap_indx], individual[i]
    # return individual,

def get_rand_mut_num():
    rand_mut = 0.03
    r = random.uniform(-rand_mut, rand_mut)
    return r

def vector_to_def(vector):
    reshaped_ind = np.array(vector).reshape(def_shape)
    return reshaped_ind

def run(out_dir, ngen):

    #Prinrt the ideal jac
    ideal_out = os.path.join(out_dir, 'ideal_jac.nrrd')
    ideal_im = sitk.GetImageFromArray(ideal)
    sitk.WriteImage(ideal_im, ideal_out)


    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()
    # pool = multiprocessing.Pool(multiprocessing.cpu_count() / 2)

    # toolbox.register("map", pool.map)

    toolbox.register("def_field", random.uniform, -0.1, 0.1)
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.def_field, n=ideal.size * 3)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    toolbox.register("evaluate", calc_fitness)
    toolbox.register("mate", cross_over)
    toolbox.register("mutate", mutate, indpb=0.1)
    toolbox.register("select", tools.selTournament, tournsize=3)

    population = toolbox.population(n=60)

    top_each_gen = []

    temp_results = os.path.join(out_dir, 'intermediate_results')
    if not os.path.exists(temp_results):
        os.mkdir(temp_results)

    intermediate_num = 50  # Spit out intermediate results every x num generations

    with open(os.path.join(out_dir, 'chart_data.txt'), 'w')as dh:

        for gen in range(ngen):
            print gen
            offspring = algorithms.varAnd(population, toolbox, cxpb=0.7, mutpb=0.1)
            fits = toolbox.map(toolbox.evaluate, offspring)
            highest = None
            for fit, ind in zip(fits, offspring):
                if not highest:
                    highest = fit
                if fit > highest:
                    highest = fit
                ind.fitness.values = fit
            top_each_gen.append(highest)
            dh.write('{}\n'.format(highest[0]))
            dh.flush()
            population = toolbox.select(offspring, k=len(population))

            if gen != 0 and gen % intermediate_num == 0:
                 top10 = tools.selBest(population, k=5)
                 write_results(temp_results, top10)


        top10 = tools.selBest(population, k=5)
        write_results(out_dir, top10)

def write_results(out, top10):
    for i, final in enumerate(top10):
        jac_result = make_jac(final)

        #write the jac
        outpath = os.path.join(out, "jac_{}.nrrd".format(i))
        sitk.WriteImage(sitk.GetImageFromArray(jac_result), outpath)

        #write the def
        outpath = os.path.join(out, "def_{}.nrrd".format(i))
        sitk.WriteImage(sitk.GetImageFromArray(vector_to_def(final)), outpath)

if __name__ == '__main__':
    import sys

    label = sys.argv[1]
    jac_value = sys.argv[2]
    out_dir = sys.argv[3]
    ngen = sys.argv[4]

    l = sitk.ReadImage(label)
    ideal = sitk.GetArrayFromImage(l)

    jac_float = float(jac_value)
    ideal[ideal == 1.0] = jac_float
    ideal[ideal == 0.0] = 1.0
    def_shape = list(ideal.shape) + [3]

    run(out_dir, int(ngen))

