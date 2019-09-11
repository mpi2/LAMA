"""
Stores default paths and has function for creating paths
"""
from os.path import join
from lama import common
from pathlib import Path
from typing import Iterator, Tuple, Dict


def specimen_iterator(reg_out_dir: Path) -> Iterator[Tuple[Path, Path]]:
    """
    Given a registration root folder/output , iterate over the line in the subfolders (could also be a single baseline folder as wel)

    Parameters
    ----------
    reg_out_dir

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

