#!/usr/bin/env python


"""
Try to implement the LM function in numpy in order to speed things up
"""


import numpy as np
import pandas as pd


#
# def run(self):
# 	data = np.vstack((self.wt_data, self.mut_data))
#
#
# def cov(x, y):
# 	X = np.column_stack([x, y]).astype(np.float32)
# 	X -= X.mean(axis=0)
# 	N = len(y)
# 	fact = N
# 	# by_hand = np.dot(X.T, X.conj()) / fact
#
# 	x1 = X.T
# 	x2 = X.conj()
# 	x_test = np.copy(x1)
# 	x_test[0][0] = 4
#
# 	cov_mat = np.matmul([x1, x1, x_test], [x2, x2, x2]) / fact
#
# 	return cov_mat
#
#
#
#
# time = df[['Time']]
# pd.DataFrame(np.linalg.pinv(time.T.dot(time)).dot(time.T).dot(df.fillna(0)),
#              ['Slope'], df.columns)

from sklearn.preprocessing import OneHotEncoder
from sklearn import preprocessing
from sklearn.linear_model import LinearRegression

mut = [5000,5100,5200,5300]
wt = [3000,3100,3200,3300]
y = wt + mut
x = [0,0,0,0,1,1,1,1]

cat_features = ['wt', 'mut']
# X is a numpy array with your features
# y is the label array
enc = OneHotEncoder(sparse=False)
X_transform = enc.fit_transform(x)

# apply your linear regression as you want
model = LinearRegression()
model.fit(X_transform, y)



# apply your linear regression as you want








