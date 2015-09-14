import numpy as np
import SimpleITK as sitk
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import rpy2.robjects as robj

rstats = importr('stats')


class AbstractStatisicalTest(object):
    def __init__(self):
        pass


class Ttest(AbstractStatisicalTest):
    def __init__(self, *args):
        super(Ttest, self).__init__(*args)


class AbstractStatsGenerator(object):
    """
    Generates the statistics. Can be all against all or each mutant against all wildtypes
    """
    def __init__(self, wt_data, mut_data, mask):
        self.wt_data = wt_data
        self.mut_data = mut_data
        self.mask = mask
        self.shape = self.set_shape()

    def run(self):
        raise NotImplementedError

    def set_shape(self):
        return self.wt_data[0].shape[0:3]


class OneAgainstMany(AbstractStatsGenerator):
    """
    Compare each mutant against all the wildtypes. Hopefully catches some variable penetrance
    """
    def __init__(self):
        pass

    def run(self):
        pass


class ManyAgainstManyStats(AbstractStatsGenerator):
    """
    Compare all the mutants against all the wild type. Generate a stats overlay
    """
    def __init__(self, *args):
        super(ManyAgainstManyStats, self).__init__(*args)
        self.stats_method_object = None

    def set_stats_method(self, stats_method_object):
        self.stats_method_object = stats_method_object

    def run(self):

        tstats, pvalues = self.stats_method_object(self.wt_data, self.mut_data)

        # print("Calculating FDR")
        qvalues = fdr(pvalues, mask_arr)

        # print 'q', min(qvalues), max(qvalues), np.mean(qvalues)
        # reshape

        filtered_t = filter_tsats(tstats, qvalues)

        t_vol = filtered_t.reshape(shape)

        t_img = sitk.GetImageFromArray(t_vol)
        outfile = os.path.join(analysis_dir, TSCORE_OUT_SUFFIX)

        sitk.WriteImage(t_img, outfile)


class AbstractFalseDiscoveryCorrection(object):
    """
    Given a set of pvalues or other statistical measure, correct based on a method defined in the subclass
    """
    def __init__(self, pvalues, mask=None):
        """
        Parameters
        ----------
        pvalues: array
            list of pvalues to correct
        mask: numpy 3D array
        """
        self.pvalues = pvalues
        self.mask = self.mask
        self.masked_pvalues = self.masked_pvalues

    def masked_pvalues(self):

        flat_mask = self.mask.flatten()
        self.pvalues[flat_mask == 0] = robj.NA_Real

    def get_qvalues(self):
        raise NotImplementedError


class BenjaminiHochberg(AbstractFalseDiscoveryCorrection):
    def __init__(self, *args):
        super(BenjaminiHochberg, self).__init__(*args)

    def get_qvalues(self):
        qvals = np.array(rstats.p_adjust(FloatVector(self.masked_pvalues), method='BH'))
        qvals[np.isnan(qvals)] = 1
        qvals[np.isneginf(qvals)] = 1
        qvals[np.isinf(qvals)] = 1
        return qvals
