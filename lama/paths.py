"""
Stores default paths and has function for creating paths
"""
from os.path import join
from lama import common
from pathlib import Path
from typing import Iterator, Tuple


def iterate_over_specimens(reg_out_dir: Path) -> Iterator[Tuple[Path, Path]]:
    """
    Given a registration output folder, iterate over the line folders (could also be a single baseline folder as wel)

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


class RegPaths(object):
    """
    Class to generate paths relative to the config file for the registration project
    If a path is not given in the config file, use a default if available
    """
    def __init__(self, config_dir_path, config):
        self.config_dir_path = config_dir_path
        self.config = config

        self.outdir = join(self.config_dir_path, 'output')

        self.default_paths = {
            'output_dir': self.outdir,
            'deformations': 'deformations',
            'jacobians': 'jacobians',
            'jacmat': 'jacobian_matrices',
            'glcms': 'glcms',
            'root_reg_dir': 'registrations',
            'inverted_labels': 'inverted_labels',
            'inverted_stats_masks': 'inverted_stats_masks',
            'organ_vol_result_csv': common.ORGAN_VOLUME_CSV_FILE
        }

    def __get_item__(self, key: str):
        """
        Just a wwrapper around get() so we can use paths['name'] syntax
        Parameters
        ----------
        key

        Returns
        -------

        """
        return self.get(key)

    def get(self, name: str) -> Path:
        """
        Get a directory path from a lama output directory

        Parameters
        ----------
        name
            The folder name

        Returns
        -------

        """
        return self.make(name, mkdir=False)

    def make(self, name, mkdir=True, parent=False):

        config_basename = self.config.get(name)

        # If the path is specified in the config use that
        if config_basename:
            if parent:
                path = join(self.outdir, parent, config_basename)
            else:
                path = join(self.outdir, config_basename)
        else:
            default_basename = self.default_paths.get(name)
            if not default_basename:
                default_basename = name
            if parent:
                path = join(self.outdir, parent, default_basename)
            else:
                path = join(self.outdir, default_basename)

        if mkdir:
            if mkdir in ('f', 'force'):
                common.mkdir_force(path)
            else:
                common.mkdir_if_not_exists(path)
        return path


