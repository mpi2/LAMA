"""
Stores default paths and has function for creating paths
"""

from pathlib import Path
from typing import Iterator, Tuple, Dict, List
from lama.elastix import REG_DIR_ORDER


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
        specimen_root = Path(specimen_root)

        if not specimen:
            self.specimen_id = specimen_root.name
        else:
            self.specimen_id = specimen

        self.reg_order = self._get_reg_order(specimen_root)
        self.outroot = specimen_root / 'output'
        self.reg_dirs = self.get_multistage_data(self.outroot / 'registrations')
        self.jacobians_dirs = self.get_multistage_data(specimen_root / 'jacobians') # Possible to have more than one
        self.deformations_dirs = self.get_multistage_data(specimen_root / 'deformations')
        self.inverted_labels_dir = self.get_multistage_data(specimen_root / 'inverted_labels')
        self.qc = specimen_root / 'output' / 'qc'
        self.qc_red_cyan_dirs = self.get_multistage_data(self.qc / 'red_cyan_overlays')
        self.red_cyan_axial = (specimen_root.name + '.png')

    def get_red_cyan_qc_images(self) -> Dict:

        qc_images = []
        for stage_dir in self.qc_red_cyan_dirs:
            if (stage_dir / 'resolution_images').is_dir():
                # Get all resolution images
                pass
            else:
                qc_images['stage']['axial'] = stage_dir / f'{self.specimen_id}_axial.png'
                qc_images['stage']['coronal'] = stage_dir / f'{self.specimen_id}_coronal.png'
                qc_images['stage']['sagittal'] = stage_dir / f'{self.specimen_id}_sagittal.png'
        return qc_images

    def get_multistage_data(self, root: Path):
        result = []
        for stage in self.reg_order:
            result.append(root / stage)
        return result

    def _get_reg_order(self, spec_root):
        """
        Text file in registrations folder that shows the orfer of registritons
        """
        order = []
        with open((spec_root / 'output' / 'registrations' / REG_DIR_ORDER), 'r') as fh:
            for line in fh:
                if line.strip():
                    order.append(line.strip())
        return order

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
                for img_path in reso_dir.iterdir():
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
        line_dir, spec_dir = self.spec_it[self.n]
        self.n += 1
        if self.n <= len(self.spec_it) -1:

            return SpecimenDataPaths(spec_dir, line_dir.name, spec_dir.name)

        else:

            raise StopIteration

    def __len__(self):
        return len(self.spec_it)


