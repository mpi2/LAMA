from statistics import mean
from numpy import std, mean, sqrt


def cohens_d(x, y):

    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2

    # For testing
    return (mean(x) - mean(y)) / sqrt(((nx - 1) * std(x, ddof=1) ** 2 + (ny - 1) * std(y, ddof=1) ** 2) / dof)