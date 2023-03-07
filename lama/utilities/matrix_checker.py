
import numpy as np


import pandas as pd
from lama.stats.heatmap import clustermap
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore
def main():

    #url = "https://raw.githubusercontent.com/dorkylever/LAMA/master/lama/tests/clustermap_data.csv"
    X = pd.read_csv(url, index_col=0)
    X.dropna(how='all', inplace=True)
    X.fillna(1, inplace=True)
    # replace missing values with 0
     # replace infinite values with 0

    print("Missing values:", X.isnull().sum().sum())
    print("Infinite values:", np.isposinf(X).sum().sum()+np.isneginf(X).sum().sum())


    #mean_data = X.apply(np.mean, axis=1)

    #std_data = X.apply(np.std, axis=1)

    #constant_rows = X.apply(lambda row: row.nunique() == 1, axis=1)

    #X = X[~constant_rows]



    #na_mean = mean_data[mean_data.isna().any()]





    #na_std = std_data[std_data.iszero().any()]

    cg = sns.clustermap(X,
                        metric="euclidean",
                        cmap=sns.diverging_palette(250, 15, l=70, s=400, sep=1, n=512, center="light",
                                                   as_cmap=True),
                        cbar_kws={'label': 'mean volume ratio'}, square=True,
                        center=1,
                        figsize=[30, len(X) * 0.3])




    # Calculate the mean and standard deviation of each variable
    #X.columns = [col.rsplit('_', 3)[-1] for col in X.columns]


    plt.tight_layout()

    plt.savefig("E:/221122_two_way/permutation_stats/rad_perm_output/two_way/easy_fig.png")
    plt.close()

    X = X.to_numpy()



    print(np.isnan(X).any() | np.isinf(X).any())
    #
    print(X)
    mu = np.mean(X, axis=0)
    sigma = np.std(X, axis=0)
    #
    #
    print(np.all(sigma))
    #
    # # Calculate the z-score matrix
    Z = (X - mu) / sigma
    print(Z)
    #
    # # Calculate the Euclidean distance matrix
    d = np.zeros((Z.shape[0], Z.shape[0]))
    for i in range(Z.shape[0]):
         for j in range(Z.shape[0]):
             d[i,j] = np.sqrt(np.sum((Z[i,:] - Z[j,:])**2))
    #
    print(d)
    print(np.isnan(d).any() | np.isinf(d).any())



if __name__ == '__main__':
    main()