"""
Here we have some classes that hld all the information about a given statistical test
- input and output paths
- References to the input and putput data
"""

from pathlib import Path


class StatsData(object):
    def __init__(self, wt_root: Path, mut_root: Path):

        self.wt_root = wt_root
        self.mut_root = mut_root

    @staticmethod
    def factory(type_):
        if type_ == 'intensity':
            return IntensityData
        if type_ == 'jacobians':
            return JacobiansData


class IntensityData(StatsData):
    def __init__(self, *args):
        super().__init__(*args)


class JacobiansData(StatsData):
    def __init__(self, *args):
        super().__init__(*args)