import multiprocessing
from deap import creator, base, tools, algorithms
import numpy as np
import random
import SimpleITK as sitk


ideal = np.ones(20*20*20).reshape((20,20,20))
ideal[2:4,2:4,2:4] = 0.7
def_shape = list(ideal.shape) + [3]



creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
pool = multiprocessing.Pool(4)

toolbox.register("map", pool.map)

toolbox.register("def_field", random.uniform, -20.0, 20.0)
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.def_field, n=ideal.size * 3)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

def calc_fitness(individual):
    jac_array = make_jac(individual)
    comp = 1 - np.sum(np.square(ideal - jac_array))
    return [comp]

def make_jac(individual):
    reshaped_ind = np.array(individual).reshape(def_shape)
    im = sitk.GetImageFromArray(reshaped_ind)
    jac = sitk.DisplacementFieldJacobianDeterminant(im)
    jac_array = sitk.GetArrayFromImage(jac)
    return jac_array


toolbox.register("evaluate", calc_fitness)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)

population = toolbox.population(n=70)

NGEN=150
for gen in range(NGEN):
    print gen
    offspring = algorithms.varAnd(population, toolbox, cxpb=0.5, mutpb=0.1)
    fits = toolbox.map(toolbox.evaluate, offspring)
    print fits
    for fit, ind in zip(fits, offspring):
        ind.fitness.values = fit
    population = toolbox.select(offspring, k=len(population))
top10 = tools.selBest(population, k=10)

for final in sorted(top10, key=lambda x: calc_fitness(x)[0]):
    jac_result = make_jac(final)
    print jac_result


