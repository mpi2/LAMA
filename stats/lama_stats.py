#!/usr/bin/env python

"""
The script that runs the statistical analysis of the LAMA pipeline. Either run standalone via the command line
or via phenodetecct.py
"""

# Hack. Relative package imports won't work if this module is run as __main__
import sys
from os.path import join, dirname, basename, abspath, splitext

sys.path.insert(0, join(dirname(__file__), '..'))

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from _phenotype_statistics import DeformationStats, IntensityStats, JacobianStats, OrganVolumeStats, AngularStats
from _stats import TTest, LinearModelR, CircularStatsTest, LinearModelPython
import common
from common import LamaDataException, Roi
import gc
import logging
from staging.get_volumes_by_stage import VolumeGetter
from stats_config_validation import validate

# Map the stats name and analysis types specified in stats.yaml to the correct class
STATS_METHODS = {
    'LM': LinearModelR,
    'LM_python': LinearModelPython,
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

DEFAULT_FORMULAS = ['genotype']  # Should add CRl as default?
DEFAULT_HEADER = ['volume_id', 'genotype', 'crl']
DEFULAT_BLUR_FWHM = 200
STAGING_PLT_NAME = 'staging.png'


def run(config_path):
    """
    The entry point to the LAMA stats script

    Parameters
    ----------
    config_path: str
        full path to the lama stats yaml config
    """

    config = validate(config_path)
    config['root_dir'] = dirname(config_path)
    config.formulas = get_formulas(config)
    wt_staging_data, mut_staging_data = get_staging_data(config.root_dir, config['wt_staging_file'], config['mut_staging_file'])

    config.wt_staging_data, config.mut_staging_data = wt_staging_data, mut_staging_data

    #  Iterate over all the stats types (eg jacobians, intensity)specified under the 'data'section of the config
    for stats_analysis_type, stats_analysis_config in config.data.iteritems():
        outdir = join(config.root_dir, stats_analysis_type)
        common.mkdir_force(outdir)
        setup_logging(outdir)
        analysis_config = stats_analysis_config
        stats_tests = analysis_config.get('tests', ['LM'])
        wt_file_list, mut_file_list, wt_basenames, mut_basenames = get_specimens(config, stats_analysis_config)
        groups_file = os.path.abspath(join(outdir, 'combined_groups.csv'))
        write_groups_file_for_r(groups_file, config, wt_basenames, mut_basenames, config.root_dir)
        staging_plot(groups_file, outdir)
        global_stats_config = setup_global_config(config) # Makes paths and sets up some defaults etc and adds back to config
        global_stats_config.groups = groups_file
        global_stats_config.wt_file_list = wt_file_list
        global_stats_config.mut_file_list = mut_file_list
        run_single_analysis(config, stats_analysis_type, outdir, stats_tests)


def setup_logging(outdir):
    """
    If there is a log file specified in the config, use that path. Otherwise log to the stats folder
    """
    logpath = join(outdir, 'stats.log')
    logging.basicConfig(filename=logpath, level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p')

    stdout_log = logging.StreamHandler(sys.stdout)
    logging.getLogger().addHandler(stdout_log)

    fileh = logging.FileHandler(logpath, 'a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fileh.setFormatter(formatter)

    log = logging.getLogger()  # root logger
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)
    log.addHandler(fileh)


def get_staging_data(root_dir, wt_staging_path, mut_staging_path):
    """
    Get the staging data for wild types and mutants


    Parameters
    ----------
    root_dir: str
        the root project directory
    wt_staging_path: str
        path to wild type staging info
    mut_staging_path
        path to mutant staging info

    Notes
    -----
        A staging file is in this format

        vol, value
        volume1.nrrd, 1.2
        volume2.nrrd, 1.3


    Returns
    -------
    tuple(dict, dict)
        wild type staging data, mutant staging data

    """
    wt_staging_data = common.csv_read_dict(join(root_dir, wt_staging_path))
    mut_staging_data = common.csv_read_dict(join(root_dir, wt_staging_path))
    return wt_staging_data, mut_staging_data


def get_specimens(global_config, stats_entry):
    """
    Get a list of wildtype and a list of mutant absolute paths to use in the analysis
    Depending on whether there is a staging 

    Parameters
    ----------
    global_config: dict
        The whole stats config object
    stats_entry: dict
        The part of the config specific to one stats run (eg. jacobians)
    plot_path: str
        path to save a plot of wt and mutant staging metric

    Returns
    -------
    pandas.DataFrame
        vol_id(minus ext), genotype (mut/wt), CRL (if available)


    TODO: Re-add the ability to specify groups files for when we have multiple effects
    """

    root_dir = global_config['root_dir']

    if stats_entry.get('wt_list'):
        wt_list_path = join(root_dir, (stats_entry['wt_list']))
        all_wt_file_list = common.get_inputs_from_file_list(wt_list_path, root_dir)

    elif stats_entry.get('wt_dir'):
        wt_data_dir = abspath(join(root_dir, stats_entry.get('wt_dir')))
        all_wt_file_list = common.GetFilePaths(wt_data_dir, ignore_folder='resolution_images')

    else:
        logging.error("A 'wt_list' or 'wt_dir' must be specified in the stats config file")
        sys.exit()

    if not all_wt_file_list:
        logging.error('Cannot find data files in {}. Check the paths in stats.yaml'.format(wt_data_dir))
        sys.exit()

    if stats_entry.get('mut_list'):
        mut_list_path = join(root_dir, (stats_entry['mut_list']))
        mut_file_list = common.get_inputs_from_file_list(mut_list_path, root_dir)

    elif stats_entry.get('mut_dir'):
        mut_data_dir = abspath(join(root_dir, stats_entry['mut_dir']))
        mut_file_list = common.GetFilePaths(mut_data_dir, ignore_folder='resolution_images')

    else:
        logging.error("A 'mut_list' or 'mut_dir' must be specified in the stats config file")
        sys.exit()

    if not mut_file_list:
        logging.error('Cannot find data files in {}. Check the paths in stats.yaml'.format(mut_data_dir))
        sys.exit()

    # Get the littermate names, if present
    littermate_file = global_config.get('littermate_controls')
    littermate_basenames = []
    if littermate_file:
        litter_mate_path = join(root_dir, littermate_file)
        littermate_basenames = common.strip_extensions(common.csv_read_lines(litter_mate_path))

    # Now we have the list of mutants and wts, if we are doing automatic staging, filter the WT list now
    wt_staging_file = global_config.get('wt_staging_file')

    # Select baselines by automatic staging unless a list of baselines is given
    if wt_staging_file and not stats_entry.get('wt_list'):
        logging.info("Choosing baselines automatically")
        mut_staging_file = global_config.get('mut_staging_file')
        if not mut_staging_file:
            logging.error("'mut_staging_file' must be specifies along with the 'wt_staging_file'")
            sys.exit(1)
        wt_file = join(root_dir, wt_staging_file)
        mut_file = join(root_dir, mut_staging_file)

        # Get the ids of volumes that are within the staging range
        mut_ids_used = common.strip_extensions([basename(x) for x in mut_file_list])
        stager = VolumeGetter(wt_file, mut_file, littermate_basenames, mut_ids_used)

        stage_filtered_wts = stager.filtered_wt_ids()

        if stage_filtered_wts is None:
            logging.error("The current staging appraoch was not able to identify enough wild type specimens")
            sys.exit(1)

        #  Keep the wt paths that were identified as being within the staging range
        wt_file_list = [x for x in all_wt_file_list
                        if basename(x).strip('seg_') in stage_filtered_wts  # filenames with extension
                        or
                        splitext(basename(x))[0].strip('seg_') in stage_filtered_wts]  # without extension

    else:  # no staging file, just use all wild types
        wt_file_list = all_wt_file_list

    wt_file_check = common.check_file_paths(wt_file_list, ret_string=True)
    if wt_file_check is not True:
        logging.error("Error: Following wild type paths for stats could not be found\n{}".format(wt_file_check))
        sys.exit(1)

    mut_file_check = common.check_file_paths(mut_file_list, ret_string=True)
    if mut_file_check is not True:
        logging.error("Error: Following mutant paths for stats could not be found\n{}".format(mut_file_check))
        sys.exit(1)

    # If we have a list of littermate basenames, remove littermates baslines from mut set and add to wildtypes
    # TODO check if littermates are in same staging range

    # remove littermates from mutant set and transfer to wt_set
    idx_to_remove = []

    if littermate_basenames:
        for lbn in littermate_basenames:
            for i in range(len(mut_file_list)):
                bn = basename(mut_file_list[i])
                bn_noext = common.strip_extension(bn)
                if lbn in (bn, bn_noext):
                    idx_to_remove.append(i)

        muts_minus_littermates = [x for i, x in enumerate(mut_file_list) if i not in idx_to_remove]
        wt_basenames = [basename(x) for x in wt_file_list]
        for idx in idx_to_remove:
            # If mut vol with same is present in wt baseline set, do not add to WT baselines.
            # RThis could happen, for instance, if the littermate controls
            # are already included in the baeline set
            if not basename(mut_file_list[idx]) in wt_basenames:
                wt_file_list.append(mut_file_list[idx])
        mut_file_list = muts_minus_littermates

    wt_basenames = [basename(x) for x in wt_file_list]  # rebuild after adding any littermates
    mut_basenames = [basename(x) for x in mut_file_list]

    if len(wt_basenames) < 1:
        logging.error("Can't find any WTs for groups file.")
        sys.exit(1)

    if len(mut_basenames) < 1:
        logging.error("Can't find any mutants for groups file.")
        sys.exit(1)

    return wt_file_list, mut_file_list, wt_basenames, mut_basenames


def staging_plot(groups_file, outdir):
    df = pd.DataFrame.from_csv(groups_file)
    sns.swarmplot(x='genotype', y='crl', data=df)
    pltpath = join(outdir, STAGING_PLT_NAME)
    plt.savefig(pltpath)


def write_groups_file_for_r(groups_file_path, config, wt_basenames, mut_basenames, root_dir):
    try:
        with open(groups_file_path, 'w') as cw:

            lm_formula = config.formulas[0]
            formula_elements = lm_formula.split(',')

            use_crl = False  # This is a bodge to get crl as a fixed effect . need to work on this
            if len(formula_elements) == 2 and config.get('wt_staging_file') and config.get('mut_staging_file'):
                use_crl = True

                wt_crls = common.csv_read_dict(join(root_dir, config.get('wt_staging_file')))
                mut_crl_file = join(root_dir, config.get('mut_staging_file'))
                mutant_crls = common.csv_read_dict(mut_crl_file)
                cw.write('volume_id,genotype,crl\n')
            else:
                cw.write('volume_id,genotype\n')

            for volname in wt_basenames:
                if use_crl:
                    cw.write('{},{},{}\n'.format(volname, 'wildtype', wt_crls[common.strip_extension(volname)]))
                else:
                    cw.write('{},{}\n'.format(volname, 'wildtype'))

            for volname in mut_basenames:
                if use_crl:
                    cw.write('{},{},{}\n'.format(volname, 'mutant', mutant_crls[common.strip_extension(volname)]))
                else:
                    cw.write('{},{}\n'.format(volname, 'mutant'))

    except (IOError, OSError) as e:
        logging.exception("Cannot open combined groups file:\n".format(groups_file_path, e.strerror))
        sys.exit(1)


def get_formulas(config):
    """
    Extract the linear model formula from the stasts config file. Just extract the independent varibale names for now

    Returns
    -------
    str: the independent variables/fixed effects
        or
    None: if no formulas can be found
    """
    parsed_formulas = []
    formulas = config.get('formulas')
    if not formulas:
        return DEFAULT_FORMULAS
    if not formulas:
        return None
    else:
        for formula_string in formulas:
            formula_elements = formula_string.split()[0::2][
                               1:]  # extract all the effects, miss out the dependent variable
            parsed_formulas.append(','.join(formula_elements))
        return parsed_formulas


def get_normalisation(config, mask_array):
    normalisation_roi = config.get('normalisation_roi')
    roi = None
    if normalisation_roi == 'mask':
        roi = mask_array
    elif isinstance(normalisation_roi, list):
        (x1, y1, z1), (x2, y2, z2) = normalisation_roi
        roi = Roi(x1=x1, x2=x2, y1=y1, y2=y2, z1=z1, z2=z2)
    elif isinstance(normalisation_roi, str):
        n = abspath(join(config.root_dir, normalisation_roi))
        if os.path.isfile(n):
            try:
                roi = common.img_path_to_array(n).ravel()
            except OSError as e:
                logging.error("Cannot read roi mask image for normalisation {}".format(n))
                sys.exit()
    return roi


def setup_global_config(config):
    """
    Build the required paths etc needed for all of the analyses and add back to the config.
    Does some checkoing for required parameters

    Parameters
    ----------
    config: addict.Dict
        The full stats config
    """
    root_dir = config['root_dir']

    mask = config.get('fixed_mask')
    if not mask:
        logging.warning('No mask specified in stats config file. A mask is required for stats analysis')
        return

    fixed_mask = config.fixed_mask = join(root_dir, config.fixed_mask)

    if not os.path.isfile(fixed_mask):
        logging.warning("Can't find mask {}. A mask is needed for the stats analysis".format(fixed_mask))
        return

    voxel_size = config.get('voxel_size')
    if not voxel_size:
        voxel_size = common.DEFAULT_VOXEL_SIZE
        logging.warn("Voxel size not set in config. Using a default of {}".format(common.DEFAULT_VOXEL_SIZE))

    config.voxel_size = float(voxel_size)

    config.project_name = config.get('project_name', '_')

    config.mask_array_flat = common.img_path_to_array(fixed_mask).ravel()

    invert_config = config.get('invert_config_file')
    if invert_config:
        config.invert_config_file = join(root_dir, invert_config)

    # Get the label maps and organ names, if used
    label_map_path = config.get('label_map_path')
    if isinstance(label_map_path, str):
        lp = abspath(join(root_dir, label_map_path))
        label_map = common.img_path_to_array(lp)
    else:
        label_map = None
    config.label_map_path = label_map

    organ_names_path = config.get('organ_names')

    if organ_names_path:
        onp = abspath(join(root_dir, organ_names_path))
        organ_names = {}
        with open(onp, 'rb') as onf:
            for i, line in enumerate(onf):
                organ_names[i + 1] = line.strip()
    else:
        organ_names = None
    config.organ_names = organ_names

    config.blur_fwhm = config.get('blur_fwhm', DEFULAT_BLUR_FWHM)

    return config


def run_single_analysis(config, analysis_name, outdir, stats_tests):
    """
    For each entry under 'data' in the config, setup the specific settings for it. Made this a separate function
    so it's easier to test
    Parameters
    ----------
    config: addict.Dict
        the main stats config
    analysis_name: str
        for example: intensity, jacobian, deformations
    """

    analysis_config = config.data[analysis_name]

    config.normalisation_roi = get_normalisation(analysis_config, config.mask_array_flat)
    gc.collect()

    logging.info('#### doing {} stats ####'.format(analysis_name))

    analysis_prefix = analysis_name.split('_')[0]
    stats_method = ANALYSIS_TYPES[analysis_prefix]

    # Change data_dir to data_paths lists
    stats_object = stats_method(outdir, analysis_prefix, config)
    for test in stats_tests:
        if test == 'LM' and not common.is_r_installed():
            logging.warning("Could not do linear model test for {}. Do you need to install R?".format(analysis_name))
            continue
        stats_object.run(STATS_METHODS[test], analysis_name)
        if config.invert_config_path:
            # Bodge: Organ volume stats is not invertable
            if analysis_prefix == 'organvolumes':
                continue
            stats_object.invert(config.invert_config_path)
    del stats_object


if __name__ == '__main__':
    # Log all uncaught exceptions
    sys.excepthook = common.excepthook_overide

    import argparse

    parser = argparse.ArgumentParser("Stats component of the phenotype detection pipeline")
    parser.add_argument('-c', '--config', dest='config', help='yaml config file contanign stats info', required=True)
    args = parser.parse_args()
    run(args.config)
