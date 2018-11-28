"""
Here we have some classes that hld all the information about a given statistical test
- input and output paths
- References to the input and putput data
"""

from pathlib import Path
import numpy as np
import datetime


class StatsData():
    def __init__(self, wt_root: Path,
                 mut_root: Path,
                 out_dir: Path,
                 mask: np.ndarray):
        """

        Parameters
        ----------
        wt_root
            Wild type folder
        mut_root
            mutant fodler
        out_dir
        mask
            3D mask
        """

        self.wt_root = wt_root
        self.mut_root = mut_root
        self.out_dir = out_dir
        self.mask = mask
        self.loader = None
        self.make_out_dir()

    @staticmethod
    def factory(type_):
        if type_ == 'intensity':
            return IntensityData
        elif type_ == 'jacobians':
            return JacobiansData
        elif type_ == 'organ_volume':
            return OrganVolumeData

    def make_out_dir(self):
        out = self.out_dir / self.type_
        out.mkdir(exits_ok=True)

    def load(self):
        """
        Load data from the associated data loader
        """
        pass


class IntensityData(StatsData):
    def __init__(self, *args):
        super().__init__(*args)


class JacobiansData(StatsData):
    def __init__(self, *args):
        super().__init__(*args)


class OrganVolumeData(StatsData):
    def __init__(self, *args):
        super().__init__(*args)
