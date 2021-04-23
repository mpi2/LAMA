import warnings
import matplotlib

# When running from a docker image, get much warnings from sklearn, pandas etc. Turn it off
warnings.filterwarnings("ignore")

from .version import __version__

#matplotlib.use('Agg')
