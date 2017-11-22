#!/usr/bin/env python


"""
Try to implement the LM function in numpy in order to speed things up
"""

import numpy as np
from scipy.stats import t as t_

sample_n = 10000

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import OneHotEncoder
n = 8
# y = np.array([0, 1,2,3, 5,6,7,8])
# y = y.reshape((8,1))
x = np.array([0,0,0,0,2,2,2,2])
# enc = OneHotEncoder()
x = x.reshape((8,1))

y = np.random.random(sample_n * 8).reshape((8, sample_n))

# y = np.column_stack((y, y))

# y[0][0] = 10

model = LinearRegression(fit_intercept=True, copy_X=True)
fit = model.fit(x, y)
pred = fit.predict(x)

# Test until I do maean on column
mean = np.mean(x)

se_slope = np.sqrt(np.sum((y-pred)**2, axis=0)/(n-2)) / np.sqrt(np.sum(((x - mean) **2), axis=0))
print('se_slope', se_slope)
coef = fit.coef_.flatten()
print('coef', coef)
t = coef/se_slope
print('t', t)
p = t_.sf(t, n-2)*2  # *2 for two sided test
print('p', p)







