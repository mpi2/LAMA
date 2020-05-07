"""
Stores default paths and has function for creating paths
"""
from os.path import join
from lama.registration_pipeline.validate_config import LamaConfig
from pathlib import Path
from typing import Iterator, Tuple, Dict, List


# TODO: Link up this code with where the folders are cerated during a LAMA run. Then when changes to folder names occur
# they are replfected in this iterator


def specimen_iterator(reg_out_dir: Path) -> Iterator[Tuple[Path, Path]]:
    """
    Given a registration output root folder , iterate over the speciemns of each line in the subfolders
    Note: lama considers the basliene as a single line.

    Parameters
    ----------
    reg_out_dir
        eg: E14.5/baselines/output (which contains one folder 'baseline')
        eg: E14.5/mutants/output (which contains multiple folder, one oer mutant line)

    Yields
    -------
    tuple:
        The path to the line directory
        The path to the specimen directory

    """

    if not reg_out_dir.is_dir():
        raise FileNotFoundError(f'Cannot find output directory {reg_out_dir}')

    for line_dir in reg_out_dir.iterdir():

        if not line_dir.is_dir():
            continue

        if str(line_dir).endswith('_'):  # Non specimen directories have _ suffix
            continue

        for specimen_dir in line_dir.iterdir():

            if not specimen_dir.is_dir():
                continue

            if str(specimen_dir).endswith('_'):
                continue

            yield line_dir, specimen_dir


class SpecimenDataPaths:
    """
    Contains paths for data output in a LAMA run. Not all data is currently included
    """
    def __init__(self, specimen_root: Path, line='', specimen=''):
        # These data are output per stage.
        self.line_id = line
        self.specimen_id = specimen
        self.reg_order = self._get_reg_order(specimen_root)
        self.reg_dirs = self.get_multistage_data(specimen_root / 'registrations')
        self.jacobians_dirs = self.get_multistage_data(specimen_root / 'jacobians') # Possible to have more than one
        self.deformations_dirs = self.get_multistage_data(specimen_root / 'deformations')
        self.inverted_labels_dir = self.get_multistage_data(specimen_root / 'inverted_labels')
        self.qc = specimen_root / 'output' / 'qc'
        self.qc_red_cyan = self.qc / 'cyan_red_overlay' / (specimen_root.name + '.png')

    def get_multistage_data(self, root: Path):
        result = []
        for stage in self.reg_order:
            result.append(root / stage)
        return result

    def _get_reg_order(self, spec_root):
        order = []
        with open((spec_root / 'output' / 'registrations' / 'reg_order.txt'), 'r') as fh:
            for line in fh:
                if line.strip():
                    order.append(line.strip())
        return order


class DataIterator:
    """
     I want to replace the specien iterator with something that return various data paths for each specimen
     Wrap this class around it for now to not break other code taht uses it.
    """
    def __init__(self, reg_out_dir):
        """
        Parameters
        ----------
        reg_out_dir
        The output of a job_runner run. TODO: what if the output of a lam_reg run?
        eg: E14.5/baselines/output (which contains one folder 'baseline')
        eg: E14.5/mutants/output (which contains multiple folder, one oer mutant line)

        Yields
        -------
        SpecimenDataPath
        """
        self.reg_out_dir = reg_out_dir
        self.spec_it = list(specimen_iterator(self.reg_out_dir))

    def __iter__(self):
        self.n = 0
        return self

    def __next__(self):
        for line_dir, spec_dir in self.spec_it:
            self.n += 1
            if self.n <= len(self.spec_it):

                return SpecimenDataPaths(spec_dir, line_dir.name, spec_dir.name)

            else:

                raise StopIteration

    def __len__(self):
        return len(self.spec_it)


