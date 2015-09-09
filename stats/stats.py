
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

    def run():
        pass


class ManyAgainstManyStats(AbstractStatsGenerator):
    """
    Compare all the mutants against all the wild type. Generate a stats overlay
    """
    def __init__(self, *args):
        super(ManyAgainstManyStats, self).__init__(*args)

    def run(self):

        tstats, pvalues = ttest(self.wt_data, self.mut_data)

        # print("Calculating FDR")
        qvalues = fdr(pvalues, mask_arr)

        # print 'q', min(qvalues), max(qvalues), np.mean(qvalues)
        # reshape

        filtered_t = filter_tsats(tstats, qvalues)

        t_vol = filtered_t.reshape(shape)

        t_img = sitk.GetImageFromArray(t_vol)
        outfile = os.path.join(analysis_dir, TSCORE_OUT_SUFFIX)

        sitk.WriteImage(t_img, outfile)