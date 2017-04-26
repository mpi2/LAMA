#!/usr/bin/env python

import yaml
from os.path import join, dirname, basename, abspath, splitext
import sys
import os
import csv

from _phenotype_statistics import DeformationStats, IntensityStats, JacobianStats, OrganVolumeStats, AngularStats
from _stats import TTest, LinearModelR, CircularStatsTest

# Hack. Relative package imports won't work if this module is run as __main__
sys.path.insert(0, join(os.path.dirname(__file__), '..'))
import common
from common import LamaDataException
import gc
import logging
import numpy as np

# Map the stats name and analysis types specified in stats.yaml to the correct class
STATS_METHODS = {
    'LM': LinearModelR,
    'tt': TTest,
    'circular_stats': CircularStatsTest
}

ANALYSIS_TYPES = {
    'intensity': IntensityStats,
    'deformations': DeformationStats,
    'jacobians': JacobianStats,
    'angles': AngularStats,
    'organvolumes': OrganVolumeStats
}

DEFAULT_FORMULAS = ['genotype']
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

    def get_groups_file_and_specimen_list(self):
        """
        The groups file is a csv that is used for the linear model analysis in R.
        Specimen lists 

        Returns
        -------
        str:
            path to groups file csv

        TODO: Re-add the ability to specify groups files for when we have multiple effects
        """

        combined_groups_file = os.path.abspath(join(self.config_dir, 'combined_groups.csv'))

        # Create default combined groups file. This is needed for running RScript for the linear model
        # Find an extry in stats.yaml to find data name
        # Get the list of ids fro the intensity directories
        for name, stats_entry in self.config['data'].iteritems():
            if stats_entry.get('wt_list'):
                wt_list_path = self.make_path(stats_entry['wt_list'])
                all_wt_file_list = common.get_inputs_from_file_list(wt_list_path, self.config_dir)
            elif stats_entry.get('wt_dir'):
                wt_data_dir = mut_data_dir = abspath(join(self.config_dir, stats_entry.get('wt_dir')))
                all_wt_file_list = common.GetFilePaths(wt_data_dir, ignore_folder='resolution_images')
            else:
                logging.error("A 'wt_list' or 'wt_dir' must be specified in the stats config file")
                sys.exit()
            if not all_wt_file_list:
                logging.error('Cannot find data files in {}. Check the paths in stats.yaml'.format(wt_data_dir))
                sys.exit()

            if stats_entry.get('mut_list'):
                mut_list_path = self.make_path(stats_entry['mut_list'])
                mut_file_list = common.get_inputs_from_file_list(mut_list_path, self.config_dir)
            elif stats_entry.get('mut_dir'):
                mut_data_dir = abspath(join(self.config_dir, stats_entry['mut_dir']))
                mut_file_list = common.GetFilePaths(mut_data_dir, ignore_folder='resolution_images')
            else:
                logging.error("A 'mut_list' or 'mut_dir' must be specified in the stats config file")
                sys.exit()
            if not mut_file_list:
                logging.error('Cannot find data files in {}. Check the paths in stats.yaml'.format(mut_data_dir))
                sys.exit()

            # Now we have the list of mutants and wts, if we are doing automatic staging filter the WT list now
            wt_staging_file = self.config.get('wt_staging_file')
            if wt_staging_file:
                mut_staging_file = self.config.get('mut_staging_file')
                if not mut_staging_file:
                    logging.error("'mut_staging_file' must be specifies along with the 'wt_staging_file'")
                    sys.exit(1)
                wt_file = self.make_path(wt_staging_file)
                mut_file = self.make_path(mut_staging_file)

                stage_filtered_wts = self.filter_filenames_by_staging(wt_file, mut_file)
                #  Keep the wt paths that were identified as being within gthe stageing range
                wt_file_list = [x for x in all_wt_file_list
                                if basename(x) in stage_filtered_wts  # filenames with extension
                                or
                                splitext(basename(x))[0] in stage_filtered_wts]  # without extension
            else:
                wt_file_list = all_wt_file_list

            # If we have a list of littermate basenames, remove littermates baslines from mut set and add to wildtypes
            # TODO check if littermates are in same staging range
            littermate_file = stats_entry.get('littermate_controls')
            if littermate_file:
                litter_mate_path = join(self.config_dir, littermate_file)
                littermate_basenames = common.csv_read_lines(litter_mate_path)
                # remove littermates from mutant set and transfer to wt_set
                idx_to_remove = []
                for lbn in littermate_basenames:
                    for i in range(len(mut_file_list)):
                        bn = basename(mut_file_list[i])
                        bn_noext = splitext(bn)[0]
                        if lbn in (bn, bn_noext):
                            idx_to_remove.append(i)

                muts_minus_littermates = [x for i, x in enumerate(mut_file_list) if i not in idx_to_remove]
                for idx in idx_to_remove:
                    wt_file_list.append(mut_file_list[idx])
                mut_file_list = muts_minus_littermates

            # TODO: remove any duplicates from the wt list. This could happen, for instance, if the littermate controls
            # also existied in the previously-created wt test set

            wt_basenames = [basename(x) for x in wt_file_list]
            mut_basenames = [basename(x) for x in mut_file_list]

            if len(wt_basenames) < 1:
                logging.error("Can't find any WTs for groups file.")
                sys.exit()

            if len(mut_basenames) < 1:
                logging.error("Can't find any mutants for groups file.")
                sys.exit()

            try:
                with open(combined_groups_file, 'w') as cw:
                    cw.write(','.join(DEFAULT_HAEDER) + '\n')
                    for volname in wt_basenames:

                        cw.write('{},{}\n'.format(volname, 'wildtype'))
                    for volname in mut_basenames:

                        cw.write('{},{}\n'.format(volname, 'mutant'))
            except (IOError, OSError):
                logging.error("Cannot open combined groups file:\n".format(combined_groups_file))
                sys.exit(1)
            break

        return combined_groups_file, wt_file_list, mut_file_list

    @staticmethod
    def filter_filenames_by_staging(wt_staging_file, mut_staging_file):
        """
        Given two staging csv files previously created by lama,
                eg:
                -----------
                wt1.nrrd,700
                wt2.nrrd,710
                wt3.nrrd,720
                wt4.nrrd,730....
                ------------
        get the list of wts that are nearest to the range of the mutants
        Parameters
        ----------
        wt_staging_file: str
            csv path with staging info (sacling factors for each id for example)
        mut_staging_file
            csv path with staging info

        Returns
        -------
        list of wild type specimen ids to use
        None if no suitable range of baselines could be found

        """
        def get_vols_from_range(df_, min_, max_):
            range_series = df_['value'].between(min_, max_)
            return list(range_series[range_series].index)

        import pandas as pd
        try:
            mut_staging_df = pd.DataFrame.from_csv(mut_staging_file)
            wt_staging_df = pd.DataFrame.from_csv(wt_staging_file)
        except Exception as e:
            logging.error('There was a problem reading in the one of the staging csv files\n{}\n{}\n\n{}'.format(
                wt_staging_file, mut_staging_file, e.message))
            sys.exit(1)

        mut_min = mut_staging_df.min().values[0]
        mut_max = mut_staging_df.max().values[0]

        wt_set = get_vols_from_range(wt_staging_df, mut_min, mut_max)

        if len(wt_set) < 8:
            # Can we get wildtypes 10% eaither size of the mut range?
            wt_set = get_vols_from_range(wt_staging_df, mut_min - (mut_min * 0.1), mut_max + (mut_max * 0.1))
            if len(wt_set) < 8:
                raise LamaDataException("Cannot find a suitable set of WT baselines using current staging files given" +
                                        "\n{}\n{} ".format(wt_staging_file, mut_staging_file))

        return [str(x) for x in wt_set] # convert to str as filename will only numbers end up as numberic types

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


    def run_stats_from_config(self):
        """
        Build the required stats classes for each data type
        """

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

        try:
            groups, wt_file_list, mut_file_list = self.get_groups_file_and_specimen_list()
        except LamaDataException as e:
            logging.info('lama enountered a problem with\n{}'.format(e))
            sys.exit()

        formulas = self.get_formulas()
        if not formulas:
            formulas = DEFAULT_FORMULAS

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
            lp = abspath(join(self.config_dir, label_map_path))
            label_map = common.img_path_to_array(lp)
        else:
            label_map = None
        organ_names_path = self.config.get('organ_names')
        if organ_names_path:
            onp = abspath(join(self.config_dir, organ_names_path))
            organ_names = {}
            with open(onp, 'rb') as onf:
                for i, line in enumerate(onf):
                    organ_names[i + 1] = line.strip()
        else:
            organ_names = None

        # loop over the types of data and do the required stats analysis
        for analysis_name, analysis_config in self.config['data'].iteritems():
            stats_tests = analysis_config.get('tests', ['LM'])

            if global_blur_fwhm:
                blur_fwhm = global_blur_fwhm
            elif analysis_config.get('blur_fwhm'):
                blur_fwhm = analysis_config.get('blur_fwhm')
            else:
                blur_fwhm = None
                logging.warning("no blur radius specified, using default")

            outdir = join(self.config_dir, analysis_name)
            normalisation_roi = analysis_config.get('normalisation_roi')
            gc.collect()

            logging.info('#### doing {} stats ####'.format(analysis_name))
            analysis_prefix = analysis_name.split('_')[0]
            stats_method = ANALYSIS_TYPES[analysis_prefix]

            # Change data_dir to data_paths lists
            stats_object = stats_method(outdir, wt_file_list, mut_file_list, project_name, mask_array_flat, groups,
                                        formulas, do_n1, (subsampled_mask, subsample), normalisation_roi, blur_fwhm,
                                        voxel_size=voxel_size,
                                        label_map=label_map, label_names=organ_names)
            for test in stats_tests:
                if test == 'LM' and not self.r_installed:
                    logging.warn("Could not do linear model test for {}. Do you need to install R?".format(analysis_name))
                    continue
                stats_object.run(STATS_METHODS[test], analysis_name)
                if invert_config_path:
                    # Bodge: Organ volume stats if not invertable
                    if analysis_prefix == 'organvolumes':
                        continue
                    stats_object.invert(invert_config_path)
            del stats_object

    def run_stats_method(self):
        pass

class StatsRun():
    """
    Contains information of a single analysis generated from the config file
    """
    def __init__(self):
        pass


if __name__ == '__main__':

    # Log all uncaught exceptions
    sys.excepthook = common.excepthook_overide

    import argparse
    parser = argparse.ArgumentParser("Stats component of the phenotype detection pipeline")
    parser.add_argument('-c', '--config', dest='config', help='yaml config file contanign stats info', required=True)
    args = parser.parse_args()
    LamaStats(args.config)


