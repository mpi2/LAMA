#!/usr/bin/env python


"""
Try to implement the LM function in numpy in order to speed things up
"""

import numpy as np
from scipy.stats import t as t_

sample_n = 10000

np.random.seed(300)

from sklearn.linear_model import LinearRegression
n = 8
# y = np.array([0, 1,2,3, 5,6,7,8])
# y = y.reshape((8,1))
genotype = np.array([0, 0, 0, 0, 1, 1, 1, 1]).reshape((8, 1))
# crl = np.array([8, 8.2, 8.1, 7.9, 10, 10.2, 10.6, 9.9]).reshape((8, 1))
# crl = np.array([0] * 8).reshape((8, 1))

# x = np.column_stack((genotype, crl))

y = np.array([1, 2, 3, 4, 5, 6, 7, 8])


y = np.column_stack((y, y))

# y[0][0] = 10

model = LinearRegression(fit_intercept=True, copy_X=True)
fit = model.fit(genotype, y)
# So we have 2 coeffiencients for each sample first is for genotype second is for crl
pred = fit.predict(genotype)


print('############ Simple linear regression')
# Simple linear regression
# Test until I do maean on column
mean_x = np.mean(genotype)

se_slope = np.sqrt(np.sum((y-pred)**2, axis=0)/(n-2)) / np.sqrt(np.sum(((genotype - mean_x) ** 2), axis=0))

print('se_slope', se_slope)
coef = fit.coef_.flatten()+1
print('coef', coef)
t = coef/se_slope
print('t', t)
p = t_.sf(t, n-2)*2  # *2 for two sided test
print('p', p)







