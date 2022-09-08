import pandas as pd
from lama.stats.permutation_stats.distributions import generate_random_combinations

def test_generate_random_combinations():
    df = pd.read_csv('/home/neil/Desktop/data.csv', index_col=0)
    generate_random_combinations(df,30)


if __name__ == '__main__':
    test_generate_random_combinations()