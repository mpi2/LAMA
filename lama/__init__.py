import warnings

# When running from a docker image, get much warnings from sklearn, pandas etc. Turn it off
warnings.filterwarnings("ignore")
