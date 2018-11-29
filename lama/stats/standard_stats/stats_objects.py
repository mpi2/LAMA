"""
Here we have some classes that hld all the information about a given statistical test
- input and output paths
- References to the input and putput data
"""

from pathlib import Path
import numpy as np
import datetime
from typing import Dict


class StatsData():
    def __init__(self, wt_root: Path, mut_root: Path, out_dir: Path, mask:np.ndarray, config: Dict):
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
        self.config = config
        self.loader = None
        self.input_data = None

    @staticmethod
    def factory(type_):
        if type_ == 'intensity':
            return IntensityData
        elif type_ == 'jacobians':
            return JacobiansData
        elif type_ == 'organ_volume':
            return OrganVolumeData

    def make_out_dir(self):
        out = self.out_dir / self.stats_type
        out.mkdir(exist_ok=True)

    def load(self):
        """
        Load data from the associated data loader
        """
        ld = self.loader(self.wt_root, self.mut_root, self.mask, self.config)
        input_data = ld.load_data()


class IntensityData(StatsData):
    def __init__(self, *args):
        super().__init__(*args)


class JacobiansData(StatsData):
    def __init__(self, *args):
        super().__init__(*args)


class OrganVolumeData(StatsData):
    def __init__(self, *args):
        super().__init__(*args)
