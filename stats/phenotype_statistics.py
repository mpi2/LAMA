from data_getters import GlcmDataGetter, DeformationDataGetter, ScalarDataGetter
from stats import ManyAgainstManyStats, OneAgainstMany

from os.path import join
import sys
import os
# Hack. Relative package imports won't work if this module is run as __main__
sys.path.insert(0, os.path.abspath('..'))
import common

class AbstractPhenotypeStatistics(object):
    """
    The base class for the statistics generators
    """
    def __init__(self, config, config_dir, outdir, mask_path):
        self.config = config
        self.config_dir = config_dir
        self.data_getter = self.set_data_getter()
        self.outdir = outdir
        self.mask = common.img_path_to_array(mask_path)

    def run(self):
        many_against_many = ManyAgainstManyStats(self.data_getter.wt_data, self.data_getter.mut_data, self.mask)
        many_against_many.run()

        # Todo: one against one analysis
        #OneAgainstMany

    def set_data_getter(self):
        raise NotImplementedError


class IntensityStats(AbstractPhenotypeStatistics):
    def __init__(self, *args):
        super(IntensityStats, self).__init__(*args)

    def set_data_getter(self):
        wt_data_dir = join(self.config_dir, self.config['wt'])
        mut_data_dir = join(self.config, self.config['mut'])
        return ScalarDataGetter(wt_data_dir, mut_data_dir)


class GlcmStats(AbstractPhenotypeStatistics):
    def __init__(self):
        pass


class JacobianStats(AbstractPhenotypeStatistics):
    def __init__(self):
        pass


class DeformationStats(AbstractPhenotypeStatistics):
    def __init__(self):
        pass





