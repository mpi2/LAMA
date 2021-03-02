"""
Stores default paths and has function for creating paths
"""

from pathlib import Path
from typing import Iterator, Tuple, Dict, List
import addict
# import lama
import os
import yaml
from lama.elastix import REG_DIR_ORDER_CFG, INVERT_CONFIG
from lama.common import cfg_load


# TODO: Link up this code with where the folders are cerated during a LAMA run. Then when changes to folder names occur
# TODO: raise error when a nonlama folder is suplied
# they are replfected in this iterator


def specimen_iterator(reg_out_dir: Path) -> Iterator[Tuple[Path, Path]]:
    """
    Given a registration output root folder , iterate over the speciemns of each line in the subfolders
    Note: lama considers the baseliene as a single line.

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

    Eaxmple folder structure showing which folder to use
    ----------------------------------------------------
    ├── baseline
    │   └── output # This folder can be used as reg_out_dir
    │       ├── baseline
    │       └── staging_info_volume.csv
    └── mutants
        └── output # This folder can be used as reg_out_dir
            ├── Ascc1
            ├── Bbox1
            ├── Copb2
            ├── Shmt2
            ├── staging_info_volume.csv
            ├── Synrg
            ├── Tm2d1
            ├── Tm2d2
            ├── Vars2
            └── Vps37d


    """
    reg_out_dir = Path(reg_out_dir)

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


class LamaSpecimenData:
    """
    Contains paths for data output in a LAMA run. Not all data is currently included just those that are used by other
    modules

    Why is setup function needed. Can it not go in the constructor?
    """
    def __init__(self, specimen_root: Path, line='', specimen='', input_dir=None):
        # These data are output per stage.
        self.line_id = line
        self.specimen_root = Path(specimen_root)

        if input_dir:
            self.input_dir = input_dir
        else:
            self.input_dir = specimen_root / 'inputs'

        if not specimen:
            self.specimen_id = specimen_root.name
        else:
            self.specimen_id = specimen

    def setup(self):
        # TODO: update this. I just moved this out of the constructor as it was failing there
        self.reg_order, self.inversion_order = self._get_reg_order(self.specimen_root)
        self.outroot = self.specimen_root / 'output'
        self.reg_dirs: Path = self.get_multistage_data(self.outroot / 'registrations')
        self.jacobians_dirs = self.get_multistage_data(self.outroot / 'jacobians')  # Possible to have more than one
        self.deformations_dirs = self.get_multistage_data(self.outroot / 'deformations')
        self.inverted_labels_dirs: Path = self.get_multistage_data(self.outroot / 'inverted_labels', self.inversion_order)
        self.qc = self.specimen_root / 'output' / 'qc'
        self.qc_red_cyan_dirs = self.qc / 'red_cyan_overlays'
        self.qc_inverted_labels = self.qc / 'inverted_label_overlay'
        self.qc_grey_dirs = self.qc / 'greyscales'
        return self

    def get_multistage_data(self, root: Path, order=None) -> List:
        """
        Given a root outdir return a list of stage-asociated data, if it exists. Else empty list

        Parameters
        ----------
        root
            specimen/output  dir

        Returns
        -------
        List of stage assocatied data directories for a specimen
        """
        result = []
        if not order:
            order = self.reg_order

        for stage in order:
            data_dir = root / stage
            if data_dir.is_dir():
                result.append(data_dir)
        return result

    def _get_reg_order(self, spec_root):
        """
        Text file in registrations folder that shows the order of registrations
        """
        reg_order = []
        inv_order = []
        with open((spec_root / 'output' / 'registrations' / REG_DIR_ORDER_CFG), 'r') as fh:
            for line in fh:
                if line.strip():
                    reg_order.append(line.strip())
        try:
            inv_order_cfg = spec_root / 'output' / 'inverted_transforms' / INVERT_CONFIG
            c = cfg_load(inv_order_cfg)
            for stage in c['inversion_order']:
                inv_order.append(stage)
        except FileNotFoundError:
            inv_order = None
        return reg_order, inv_order


    def registration_imgs(self) -> Iterator[Tuple[str, Path]]:
        """
        Generate the series of images in the order they were created in the pipeline

        Yields
        -------
        Path to img
        """
        for stage_dir in self.reg_dirs:
            stage_dir / self.specimen_id

            # If we have individual resolution imgs, output these
            reso_dir = stage_dir / self.specimen_id / 'resolution_images'
            if reso_dir.is_dir():
                for img_path in sorted(reso_dir.iterdir()):
                    yield  stage_dir.name, img_path
            # If no resolution images, just output the final registrated image for that stage
            else:
                yield stage_dir.name, stage_dir / self.specimen_id / (self.specimen_id + '.nrrd')


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

        if self.n <= len(self.spec_it) - 1:
            line_dir, spec_dir = self.spec_it[self.n]
            self.n += 1
            return LamaSpecimenData(spec_dir, line_dir.name, spec_dir.name)

        else:

            raise StopIteration

    def __len__(self):
        return len(self.spec_it)


def get_specimen_dirs(root: Path, depth=4) -> List[LamaSpecimenData]:
    # Identify all lama directoris by getting the log files
    # lama_logs = root.rglob('**/LAMA.log')

    specimen_dirs = []

    for log in [x for x in walk(root, depth) if x.name == 'LAMA.log']:
        root = log.parent
        # Take a guess at the line, probably the name of the spec dir parent
        line = root.parent.name
        s = LamaSpecimenData(log.parent, line=line)
        s.setup()
        specimen_dirs.append(s)

    return specimen_dirs


def walk(root: Path, depth=None):
    """
    Do a recursice walk for files up to a maximum depth
    """
    root = Path(root)

    if depth and depth == 1:
        for filename in root.iterdir():
            yield filename
    else:
        root_depth = len(root.parts)

        for dirpath, dirnames, filenames in os.walk(root):
            current_depth = len(Path(dirpath).parts)

            if current_depth - root_depth > depth:
                dirnames[:] = [] # Stop looking in subdirs
                continue

            for filename in filenames:
                yield Path(dirpath) / filename