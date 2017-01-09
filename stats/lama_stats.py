#!/usr/bin/env python

import yaml
from os.path import join, dirname, basename
import sys
import os
import csv

from _phenotype_statistics import DeformationStats, GlcmStats, IntensityStats, JacobianStats, OrganVolumeStats, AngularStats
from _stats import TTest, LinearModelR, CircularStatsTest

# Hack. Relative package imports won't work if this module is run as __main__
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common
import gc
import logging

# Map the stats name and analysis types specified in stats.yaml to the correct class
STATS_METHODS = {
    'LM': LinearModelR,
    'ks': TTest,
    'circular_stats': CircularStatsTest
}

ANALYSIS_TYPES = {
    'intensity': IntensityStats,
    'deformations': DeformationStats,
    'jacobians': JacobianStats,
    'glcm': GlcmStats,
    'angles': AngularStats,
    'organvolumes': OrganVolumeStats
}

DEFAULT_FORMULA = 'data ~ genotype'
DEFAULT_HAEDER = ['volume_id', 'genotype']


class LamaStats(object):
    """
    Takes a stats.yaml config file and creates appropriate PhenotypeStatitics subclasses based on which analysis is to
    be performed
    """
    def __init__(self, config_path):

        self.config_dir = dirname(config_path)
        self.config_path = config_path
        self.config = self.get_config(config_path)
        self.setup_logging()
        self.r_installed = True
        self.run_stats_from_config()

    def setup_logging(self):
        """
        If there is a log file specified in the config, use that path. Otherwise log to the stats folder
        """
        logpath = self.config.get('log')
        if not logpath:
            logpath = join(self.config_dir, 'stats.log')

        common.init_logging(logpath)
        logging.info('##### Stats started #####')
        logging.info(common.git_log())

    def make_path(self, path):
        """
        All paths are relative to the config file dir.
        Return relative paths to the config dir
        """
        return join(self.config_dir, path)

    def get_config(self, config_path):
        """
        Get the config and check for paths
        """
        try:
            with open(config_path) as fh:
                try:
                    config = yaml.load(fh)
                except Exception as e:  # Couldn't catch scanner error from Yaml
                    logging.error('Error reading stats yaml file\n\n{}'.format(e))
                    sys.exit()
        except IOError:
            logging.error("cannot find or open stats config file: {}".format(config_path))
            sys.exit(1)
        try:
            data = config['data']
        except KeyError:
            logging.error("stats config file needs a 'data' entry. Are you usin gthe correct config file?")

        return config

    def get_groups(self, wt_subset, mut_subset):
        """
        Combine group info from both the wildtype and mutants. Write out a combined groups csv file.
        If wt file basenames not specified in wt_subset_file, remove them

        Returns
        -------
        None: if no file can be found
        Dict: if file can be found {volume_id: {groupname: grouptype, ...}...}
        """
        wt_groups = mut_groups = None

        wt_g = self.config.get('wt_groups')
        mut_g = self.config.get('mut_groups')
        if all((wt_g, mut_g)):

            wt_groups = join(self.config_dir, self.config['wt_groups'])
            mut_groups = join(self.config_dir, self.config['mut_groups'])

            if not all((os.path.isfile(wt_groups), os.path.isfile(mut_groups))):
                wt_groups = mut_groups = None
                logging.warn("Can't find the wild type groups file, using default")

        combined_groups_file = os.path.abspath(join(self.config_dir, 'combined_groups.csv'))

        if all((wt_groups, mut_groups)):  # Generate the combined groups file from the given wt and mut files
            try:
                with open(wt_groups, 'r') as wr, open(mut_groups, 'r') as mr, open(combined_groups_file, 'w') as cw:
                    wt_reader = csv.reader(wr)
                    first = True
                    for row in wt_reader:
                        if first:
                            header = row
                            first = False
                            cw.write(','.join(header) + '\n')
                        else:
                            cw.write(','.join(row) + '\n')

                    reader_mut = csv.reader(mr)
                    first = True
                    for row in reader_mut:
                        if first:
                            header_mut = row
                            if header != header_mut:
                                logging.warn("The header for mutant and wildtype group files is not identical. Creating default groups file")
                                return None
                            first = False
                        else:
                            cw.write(','.join(row) + '\n')
            except (IOError, OSError):
                logging.error("Cannot open one or more of the groups files:\n{}\n".format(wt_groups, mut_groups))
                sys.exit(1)
        else:  # Create default combined groups file. This is needed for running RScript for the linear model
            # Find an extry in stats.yaml to find data name
            # Get the list of ids fro the intensity directories
            for name, stats_entry in self.config['data'].iteritems():
                wt_data_dir = join(self.config_dir, stats_entry['wt'])
                mut_data_dir = join(self.config_dir, stats_entry['mut'])
                wt_file_list = common.GetFilePaths(wt_data_dir, ignore_folder='resolution_images')
                mut_file_list = common.GetFilePaths(mut_data_dir, ignore_folder='resolution_images')
                if not wt_file_list:
                    logging.error('Cannot find data files in {}. Check the paths in stats.yaml'.format(wt_data_dir))
                    sys.exit(1)
                if not mut_file_list:
                    logging.error('Cannot find data files in {}. Check the paths in stats.yaml'.format(mut_data_dir))
                    sys.exit(1)
                wt_basenames = [basename(x) for x in common.GetFilePaths(wt_data_dir, ignore_folder='resolution_images')]
                mut_basenames = [basename(x) for x in common.GetFilePaths(mut_data_dir, ignore_folder='resolution_images')]

                try:
                    with open(combined_groups_file, 'w') as cw:
                        cw.write(','.join(DEFAULT_HAEDER) + '\n')
                        for volname in wt_basenames:
                            if wt_subset:
                                if os.path.splitext(volname)[0] in wt_subset:
                                    cw.write('{},{}\n'.format(volname, 'wildtype'))
                            else:
                                cw.write('{},{}\n'.format(volname, 'wildtype'))
                        for volname in mut_basenames:
                            if mut_subset:
                                if os.path.splitext(volname)[0] in mut_subset:
                                    cw.write('{},{}\n'.format(volname, 'mutant'))
                            else:
                                cw.write('{},{}\n'.format(volname, 'mutant'))
                except (IOError, OSError):
                    logging.error("Cannot open combined groups file:\n".format(combined_groups_file))
                    sys.exit(1)
                break

        return combined_groups_file

    def get_formulas(self):
        """
        Extract the linear/mixed model from the stasts config file. Just extract the independent varibale names for now

        Returns
        -------
        str: the independent variables/fixed effects
            or
        None: if no formulas can be found
        """
        parsed_formulas = []
        formulas = self.config.get('formulas')
        if not formulas:
            return None
        else:
            for formula_string in formulas:
                formula_elements = formula_string.split()[0::2][1:]  # extract all the effects, miss out the dependent variable
                parsed_formulas.append(','.join(formula_elements))
            return parsed_formulas

    def get_subset_list(self, subset_file):
        """
        Trim the files found in the wildtype input directory to thise in the optional subset list file
        """
        wt_vol_ids_to_use = []
        try:
            with open(subset_file, 'r') as reader:
                for line in reader:
                    vol_name = line.strip()
                    wt_vol_ids_to_use.append(vol_name)
            return wt_vol_ids_to_use
        except (OSError, IOError):
            logging.error("Cannot find specimen subset file: {}".format(subset_file))
            sys.exit(1)

    def get_subset_ids(self):
        """
        Get the subset list of vol ids to do stats with

        Returns
        -------
        tuple
            None if no subset file specified
            list of ids
        """
        wt_subset_file = self.config.get('wt_subset_file')
        mut_subset_file = self.config.get('mut_subset_file')
        wt_subset_ids = mut_subset_ids = None
        if wt_subset_file:
            wt_subset_file = join(self.config_dir, wt_subset_file)
            wt_subset_ids = self.get_subset_list(wt_subset_file)
            if len(wt_subset_ids) < 1:
                wt_subset_ids = None
        if mut_subset_file:
            mut_subset_file = join(self.config_dir, mut_subset_file)
            mut_subset_ids = self.get_subset_list(mut_subset_file)
            if len(mut_subset_ids) < 1:
                mut_subset_ids = None

        return wt_subset_ids, mut_subset_ids

    def run_stats_from_config(self):
        """
        Build the required stats classes for each data type
        """

        wt_subset_ids, mut_subset_ids = self.get_subset_ids()

        if self.config.get('blur_fwhm'):
            global_blur_fwhm = self.config.get('blur_fwhm')  # Apply this width of guassain to all the data sets
        else:
            global_blur_fwhm = None

        mask = self.config.get('fixed_mask')
        if not mask:
            logging.warn('No mask specified in stats config file. A mask is required for stats analysis')
            return
        fixed_mask = self.make_path(self.config.get('fixed_mask'))
        if not os.path.isfile(fixed_mask):
            logging.warn("Can't find mask {}. A mask is needed for the stats analysis".format(fixed_mask))
            return

        voxel_size = self.config.get('voxel_size')
        if not voxel_size:
            voxel_size = 28.0
            logging.warn("Voxel size not set in config. Using a default of 28")
        voxel_size = float(voxel_size)

        groups = self.get_groups(wt_subset_ids, mut_subset_ids)

        formulas = self.get_formulas()
        if not formulas:
            formulas = ['data ~ genotype']  # Default formula for Linear model

        project_name = self.config.get('project_name')
        if not project_name:
            project_name = '_'
        do_n1 = self.config.get('n1')

        mask_array = common.img_path_to_array(fixed_mask)
        mask_array_flat = mask_array.ravel()

        subsample = self.config.get('subsample')
        if subsample:
            try:
                int(subsample)
            except ValueError:
                logging.warn("Subsampling of stats not carried out. subsample: subsample factor(integer) not corretly specified in config file")
                subsample = False
                subsampled_mask = None
            else:
                subsampled_mask = common.subsample(mask_array, subsample, mask=True).ravel()
        else:
            subsampled_mask = None

        invert_config = self.config.get('invert_config_file')
        if invert_config:
            invert_config_path = self.make_path(invert_config)
        else:
            invert_config_path = None

        # Get the label maps and organ names, if used
        label_map_path = self.config.get('label_map_path')
        if label_map_path:
            lp = join(self.config_dir, label_map_path)
            label_map = common.img_path_to_array(lp)
        else:
            label_map = None
        organ_names_path = self.config.get('organ_names')
        if organ_names_path:
            onp = join(self.config_dir, organ_names_path)
            organ_names = {}
            with open(onp, 'rb') as onf:
                for i, line in enumerate(onf):
                    organ_names[i + 1] = line.strip()
        else:
            organ_names = None

        # loop over the types of data and do the required stats analysis
        for analysis_name, analysis_config in self.config['data'].iteritems():
            stats_tests = analysis_config['tests']
            if global_blur_fwhm:
                blur_fwhm = global_blur_fwhm
            elif analysis_config.get('blur_fwhm'):
                blur_fwhm = analysis_config.get('blur_fwhm')
            else:
                blur_fwhm = None
                logging.warn("no blur radius specified, using default")

            mut_data_dir = self.make_path(analysis_config['mut'])
            wt_data_dir = self.make_path(analysis_config['wt'])
            outdir = join(self.config_dir, analysis_name)
            normalisation_roi = analysis_config.get('normalisation_roi')
            gc.collect()

            logging.info('#### doing {} stats ####'.format(analysis_name))
            analysis_prefix = analysis_name.split('_')[0]
            stats_method = ANALYSIS_TYPES[analysis_prefix]
            stats_object = stats_method(outdir, wt_data_dir, mut_data_dir, project_name, mask_array_flat, groups, formulas, do_n1,
                          (subsampled_mask, subsample), normalisation_roi, blur_fwhm,
                                        voxel_size=voxel_size, wt_subset=wt_subset_ids, mut_subset=mut_subset_ids,
                                        label_map=label_map, label_names=organ_names)
            for test in stats_tests:
                if test == 'LM' and not self.r_installed:
                    logging.warn("Could not do linear model test for {}. Do you need to install R?".format(analysis_name))
                    continue
                stats_object.run(STATS_METHODS[test], analysis_name)
                if invert_config_path:
                    stats_object.invert(invert_config_path)
            del stats_object

    def run_stats_method(self):
        pass

if __name__ == '__main__':

    # Log all uncaught exceptions
    sys.excepthook = common.excepthook_overide

    import argparse
    parser = argparse.ArgumentParser("Stats component of the phenotype detection pipeline")
    parser.add_argument('-c', '--config', dest='config', help='yaml config file contanign stats info', required=True)
    args = parser.parse_args()
    LamaStats(args.config)

