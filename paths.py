"""
Stores default paths and has function for creating paths
"""
from os.path import join
import common



class RegPaths(object):
    """
    Class to generate paths relative to the config file for the registration project
    If a path is not given in the config file, use a default if available
    """
    def __init__(self, config_dir_path, config):
        self.config_dir_path = config_dir_path
        self.config = config

        if config.get('output_dir'):
            self.default_outdir = join(self.config_dir_path, config.get('output_dir'))
        else:
            self.default_outdir = join(self.config_dir_path, 'output')

        self.default_paths = {
            'output_dir': self.default_outdir,
            'deformations': 'deformations',
            'jacobians': 'jacobians',
            'jacmat': 'jacobian_matrices',
            'glcms': 'glcms',
            'root_reg_dir': 'registrations',
            'inverted_labels': 'inverted_labels',
            'inverted_masks': 'inveretd_masks'
        }

    def get(self, name):
        return self.make(name, mkdir=False)

    def make(self, name, mkdir=True, parent=False):

        config_basename = self.config.get(name)

        if config_basename:
            if parent:
                path = join(self.default_outdir, parent, config_basename)
            else:
                path = join(self.default_outdir, config_basename)
        else:
            default_basename = self.default_paths.get(name)
            if not default_basename:
                default_basename = name
            if parent:
                path = join(self.default_outdir, parent, default_basename)
            else:
                path = join(self.default_outdir, default_basename)

        if mkdir:
            if mkdir in ('f', 'force'):
                common.mkdir_force(path)
            else:
                common.mkdir_if_not_exists(path)
        return path


