#!/usr/bin/env python3

"""
The script that runs the statistical analysis of the LAMA pipeline. Either run standalone via the command line
or via jobrunner.py

This is the series of scripts that are run

    run_lama_stats.py
        This script get's the data ready from statistical analysis. Baseline comparison volumes are selected. For each
        stats entry in the config yaml file the correct stats runner class is called in analysis_types.py

    analysis_types.py
        Each class in here runs a different type of data such as jacobians intensity etc.

    statistical_tests.py
        The actual tests are run here


Notes
-----

auto_staging is False by default (All baselines will be used)

"""


import sys
import os
from lama import common
common.disable_warnings_in_docker()
from typing import Union
import shutil
from os.path import join, dirname, abspath, basename, splitext
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import gc
from logzero import logger as logging

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from lama.stats.analysis_types import DeformationStats, IntensityStats, JacobianStats, OrganVolumeStats, GlcmStats
from lama.stats.statistical_tests import LinearModelR, LinearModelNumpy
from lama.common import LamaDataException, Roi
from lama.staging.baseline_selection import BaselineSelector
from lama.stats.stats_config_validation import validate


# Map the stats name and analysis types specified in stats.yaml to the correct class
STATS_METHODS = {
    'LM': LinearModelR,
    'LM_python': LinearModelNumpy
}

ANALYSIS_TYPES = {
    'intensity': IntensityStats,
    'deformations': DeformationStats,
    'jacobians': JacobianStats,
    'organvolumes': OrganVolumeStats,
    'glcm': GlcmStats
}

DEFAULT_FORMULAS = ['genotype,crl']  # Should add CRl as default?
DEFAULT_HEADER = ['volume_id', 'genotype', 'crl']
DEFULAT_BLUR_FWHM: int = 100
STAGING_PLT_NAME: str = 'staging.png'


old_style_bodge = False



class TestData(object):
    """

    """
    def __init__(self):
        pass


def run(config_path):
    """
    The entry point to the LAMA stats script

    Parameters
    ----------
    config_path: str
        full path to the lama stats yaml config
    """

    setup_logging(Path(config_path).parent)  # Log to main stats log in the root directory

    try:
        config = validate(config_path)
    except ValueError as e:
        logging.exception("Problem reading the stats config file. See the stats.log file")
        sys.exit()
    except IOError as e:
        logging.exception(f"Problem with some paths. See the stats.log file")
        raise

    config['root_dir'] = dirname(config_path)  # default is to have root dir the same as stats config dir
    root_dir = config['root_dir']
    out_dir = config.get('output_dir')
    if out_dir:
        config.outdir = join(config.root_dir, out_dir)
    else:
        config.outdir = config.root_dir

    def get_abs_path_from_config(filename):
        if config.get(filename):
            return join(root_dir, config.get(filename))
        else:
            return None

    # 160818. I should do this for all the paths in one place
    if config.get('line_calibrated_p_values'):
        config.line_calibrated_p_values = get_abs_path_from_config('line_calibrated_p_values')
    else:
        config.line_calibrated_p_values = None
    if config.get('specimen_calibrated_p_values'):
        config.specimen_calibrated_p_values = get_abs_path_from_config('specimen_calibrated_p_values')
    else:
        config.specimen_calibrated_p_values = None

    config.formulas = get_formulas(config)

    # if config.get('use_auto_staging') and not all([config.get('wt_staging_file'), config.get('mut_staging_file')]):
    #     raise ValueError(""""Cannot do auto staging as there are no
    #     'mut_staging_file' and 'wt_staging_file' specified in the stast config""")

    if config.get('use_auto_staging') is not False:
        auto_staging = True

    else:
        auto_staging = False

    wt_staging_data = get_staging_data(Path(config.root_dir) / config['wt_root'])
    mut_staging_data = get_staging_data(Path(config.root_dir) / config['mut_root'])

    config.label_map, config.label_names = \
        get_labels_and_names(Path(config.root_dir), Path(config.get('label_map')), Path(config.get('label_names')))

    #  Iterate over all the stats types (eg jacobians, intensity) specified under the 'data'section of the config
    for stats_analysis_type, stats_analysis_config in config.data.items():
        try:
            outdir = join(config.outdir, stats_analysis_type)
            common.mkdir_force(outdir)

            # Start logging in the stats run folder
            setup_logging(outdir)

            #copy the stats config into the stats output directory just for recording what has been done
            shutil.copy(config_path, outdir) # Why are there out_dir and outdir valiables in the same scope?

            logging.info('Doing stats for {}'.format(stats_analysis_type))
            analysis_config = stats_analysis_config
            stats_tests = analysis_config.get('tests', ['LM'])

            all_wt_paths = get_file_paths(root_dir, config['wt_root'], stats_analysis_config['wt'])
            if not all_wt_paths:
                raise LamaDataException("Cannot find the wild type file paths using wt:{}".format(stats_analysis_config['wt']))

            all_mut_paths = get_file_paths(root_dir,config['mut_root'], stats_analysis_config['mut'])
            if not all_mut_paths:
                raise LamaDataException("Cannot find the mutant file paths using mut:{}".format(stats_analysis_config['mut']))

            littermates = config.get('littermate_controls')
            if littermates:
                littermates = common.strip_img_extensions(littermates)

            littermate_pattern = config.get('littermate_pattern')

            # Subset of speciemns to use. If None, use all
            mutant_ids = config.get('mutant_ids')
            wildtype_ids = config.get('wildtype_ids')

            filtered_wts, filtered_muts = get_filtered_paths(all_wt_paths,
                                                             all_mut_paths,
                                                             outdir,
                                                             mutant_ids,
                                                             wildtype_ids,
                                                             littermates,
                                                             littermate_pattern,
                                                             wt_staging_data,
                                                             mut_staging_data,
                                                             auto_staging)

            wt_basenames = [basename(x) for x in filtered_wts]
            mut_basenames = [basename(x) for x in filtered_muts]

            groups_file = os.path.abspath(join(outdir, 'combined_groups.csv'))
            write_groups_file_for_r(groups_file, wt_basenames, mut_basenames, wt_staging_data, mut_staging_data)

            if auto_staging:
                staging_plot(groups_file, outdir)

            # TODO: what if there are no label map or names?

            # Make paths and sets up some defaults etc and add back to config
            global_stats_config = setup_global_config(config) # I've forgot what global_stats_config does
            global_stats_config.groups = groups_file
            global_stats_config.wt_file_list = filtered_wts
            global_stats_config.mut_file_list = filtered_muts
            run_single_analysis(config, stats_analysis_type, outdir, stats_tests)
        except (ValueError, IOError) as e:  # Catch the error here so we can move on to next analysis if need be
            print(('stats failed for {}. See log file'.format(stats_analysis_type)))
            logging.exception('Stats failed for {}'.format(stats_analysis_type))
            logging.exception(e)
            raise


def setup_logging(outdir):

    logpath = Path(outdir) / 'stats.log'
    common.init_logging(logpath)


def get_staging_data(root_dir) -> pd.DataFrame:
    """
    Each specimen output folder will contain a staging file.
    This csv file has two columns:
        vol: the volume name minus extension
        value: the staging value. For example affine scaling factor
               which apprximates crown to rump relative to the target

    collate all these value into one dataframe.


    Parameters
    ----------
    root_dir: str
        the root project directory containing eith mutant or baselines registrations

    Returns
    -------
    pd.Dataframe
        Concatenated staging info, one row per specimen
        colums:
            vol: vol_id (str)
            value: staging metric (float)
    """

    staging_files = list(Path(root_dir).glob(f'**/{common.STAGING_INFO_FILENAME}'))
    dfs = []
    for file_ in staging_files:
        df = pd.read_csv(file_, index_col=0)
        if 'vol' not in df:
            df['vol'] = df.index
            df.columns = ['value', 'vol']

        dfs.append(df)

        if old_style_bodge:
            break #260918 remove this is just to get it working for old datastructure

    staging_df = pd.concat(dfs)

    return staging_df


def get_file_paths(project_root: str, input_root: str, path:str) -> list:
    """
    ----------
    project_root:
        path to the stats project directory (where the stats.yaml config resides)

    input_root:
        The root directory where the input data is.

        This could be:
            The folder containing all the baselines individual runs
            The folder containing a single mutant run

    path:
        The realtive path where the data of interest can be found
        eg: output/jacobians/deformable_stage
    organ_volume_csv
        If true search for partial path and if found join to project root.
        No need to then search for images in subdirectories


    Returns
    -------
    A list of absolute paths to volumes

    """

    input_dir = Path(project_root) / input_root

    dir_src = f'**{os.path.sep}{Path(path)}'
    dirs = list(input_dir.glob(str(dir_src)))

    files = []

    for d in dirs:
       files.extend(common.get_file_paths(d, ignore_folder='resolution_images'))


    # if isdir(data_dir):
    #     file_list = common.get_file_paths(data_dir,  ignore_folder='resolution_images')
    #     if not file_list:
    #         return None
    # else:
    #     file_list = common.get_inputs_from_file_list(data_dir, root_dir)
    #     if not file_list:
    #         return None

    # elif dir_or_path_file.__iter__():
    #     if len(dir_or_path_file) < 2:
    #         return None
    #     root = dir_or_path_file[0]
    #     dir_name = dir_or_path_file[1]




    return files


def filter_specimens_by_id(specimens: list, ids_to_include:Union[None, list]):
    """

    Parameters
    ----------
    specimens:
        the paths of all the specimens
    ids_to_include:
        list: (filenames of speciemns to include. can be with or without file extension)
        None: No predfined list of specimens to use

    Returns
    -------
    list
        subset of speciemns whose ID is in ids_to_include

    """
    if not ids_to_include:
        return specimens
    # in case ids are only digits, convert to string
    ids_to_include = [str(common.specimen_id_from_file_path(x)) for x in ids_to_include]
    to_use = [x for x in specimens if common.specimen_id_from_file_path(x) in ids_to_include]

    all_specimen_ids = common.specimen_ids_from_paths(specimens)  # This strips file extension
    ids_to_use_not_in_specimens = \
        set(ids_to_include).difference(set(all_specimen_ids))

    if len(ids_to_use_not_in_specimens) > 0:
        raise ValueError('\n\n{}\n is/are ids listed in the config to include in analysis, '
                         'but was not found in the following specimen list\n{}'.format(
                          '\n'.join(list(ids_to_use_not_in_specimens)), '\n'.join(specimens)))
    return to_use


def get_filtered_paths(wildtypes,
                       mutants,
                       out_dir,
                       mutant_ids_to_include=None,
                       wildtype_ids_to_include=None,
                       littermate_controls=None,
                       littermate_pattern=None,
                       wildtype_staging=None,
                       mutant_staging=None,
                       auto_staging=False):
    """

    Using various critea, create a final list of wild types and mutants to use in the analysis

    Parameters
    ----------
    wildtypes: list
        all the wild types in the given folder
    mutants: list
        all the mutants in the given folder
    out_dir: str
        path to the stats outdir. Used for writeinting a warning file is required
    mutant_ids_to_include: list
        mutants to include in the analysis. If len > 0, only mutantd in this list will be used
    wildtype_ids_to_include: list
        wild types to include in the analysis. If len > 0, only mutantd in this list will be used
    littermate_controls: list
        ids of volumes that are in the mutant set, but are actually littermate controls
        add these to the wild types
    littermate_pattern: str
        any mutant with this string with the filename (eg: _WT_) will be added to the wild types
    wildtype_staging: pd.Dataframe
        columns:
            vol, value
    mutant_staging: pd.Dataframe
         columns:
            vol, value

    Returns
    -------
    wild type file list
    mutant file list

    """
    mutants = filter_specimens_by_id(mutants, mutant_ids_to_include)
    wildtypes = filter_specimens_by_id(wildtypes, wildtype_ids_to_include)

    # Get the littermate names, if present
    littermate_basenames = []
    if littermate_pattern:
        for mut_file in mutants:
            if littermate_pattern.lower() in mut_file.lower():
                littermate_basenames.append(common.strip_img_extension(mut_file))

    # Add littermate controls to the WTs from the mutants
    if isinstance(littermate_controls, list):
        littermate_basenames.extend(common.strip_img_extensions(littermate_controls))

    # Select baselines by automatic staging unless a list of baselines is given
    if auto_staging:

        # Get the ids of volumes that are within the staging range
        mutant_basenames = common.strip_img_extensions([basename(x) for x in mutants])
        stager = BaselineSelector(wildtype_staging, mutant_staging, littermate_basenames, mutant_basenames)

        ignore_constraint = False if auto_staging else True
        stage_filtered_wts = stager.filtered_wt_ids(ignore_constraint=ignore_constraint)
        #littermate_ids_to_add_to_baselines = stager.littermates_to_include()
        littermate_ids_to_add_to_baselines = []  #100718 do not use littermates as it breaks organvolume stats

        excluded_mutants = stager.excluded_mutants
        if excluded_mutants:
            temp = []
            for m in mutants:
                if common.specimen_id_from_file_path(m) in excluded_mutants:
                    logging.info('Excluded {} due to size deviation from baselines'.format(common.specimen_id_from_file_path(m)))
                else:
                    temp.append(m)
            mutants = temp

        if stage_filtered_wts is None:
            logging.error("The current staging appraoch was not able to identify enough wild type specimens")
            logging.warn("******\nSelecting baselines that are nearest in developmental proxy\n"
                         "This may increase chances of spurious results\n******")

            open(join(out_dir, 'baseline_selection_WARNING_see_log'), 'a').close()

            stage_filtered_wts = stager.filtered_wt_ids(ignore_constraint=True)
        if stage_filtered_wts is None:
            raise LamaDataException("No baselines could be found")

        #  Keep the wt paths that were identified as being within the staging range
        wt_file_list = [x for x in wildtypes
                        if basename(x).strip('seg_') in stage_filtered_wts  # filenames with extension
                        or
                        splitext(basename(x))[0].strip('seg_') in stage_filtered_wts]  # without extension

    else:  # no staging file, just use all wild types
        wt_file_list = wildtypes
        littermate_ids_to_add_to_baselines = [] # Not using littermate WTs for now

    wt_file_check = common.check_file_paths(wt_file_list, ret_string=True)
    if wt_file_check is not True:
        logging.error("Error: Following wild type paths for stats could not be found\n{}".format(wt_file_check))
        raise LamaDataException("Error: Following wild type paths for stats could not be found\n{}".format(wt_file_check))

    mut_file_check = common.check_file_paths(mutants, ret_string=True)
    if mut_file_check is not True:
        logging.error("Error: Following mutant paths for stats could not be found\n{}".format(mut_file_check))
        raise LamaDataException("Error: Following mutant paths for stats could not be found\n{}".format(mut_file_check))

    # If we have a list of littermate basenames, remove littermates baslines from mut set and add to wildtypes
    # TODO check if littermates are in same staging range

    # remove littermates from mutant set and add to wts
    if littermate_basenames:
        for lbn in littermate_basenames:
            for mut in mutants:
                if common.strip_img_extension(basename(lbn)) == common.strip_img_extension(basename(mut)):
                    mutants.remove(mut)  # Remove liitermates from the baselines
                    if littermate_ids_to_add_to_baselines and common.specimen_id_from_file_path(mut) in littermate_ids_to_add_to_baselines:  # If within mutat CRL range add to baseline set
                        wt_file_list.append(mut)

    # If mut vol with same name is present in wt baseline set, do not add to WT baselines.
    # This could happen, for instance, if the littermate controls are already included in the baseline set
    wt_file_list = list(set(wt_file_list))

    if len(mutants) < 1:
        logging.error("Can't find any WTs for groups file.")  # What does this mean
        raise LamaDataException("Can't find any WTs for groups file.")

    if len(wt_file_list) < 1:
        logging.error("Can't find any mutants for groups file.")
        raise LamaDataException("Can't find any mutants for groups file.")

    return wt_file_list, mutants


def staging_plot(groups_file, outdir):
    df = pd.read_csv(groups_file)
    df.crl = df.crl.astype(float)
    sns.swarmplot(x='genotype', y='crl', data=df)
    pltpath = join(outdir, STAGING_PLT_NAME)
    plt.savefig(pltpath)
    plt.close()


def write_groups_file_for_r(groups_file_path:str, wt_basenames:list, mut_basenames:list, wt_staging=None, mut_staging=None):
    """
    Write out a csv that is used by the R script to run the linear model.

    The outpuit file should loook something like this

        volume_id,genotype,crl
        test.nrrd,wildtype,0.97
        test1.nrrd,mutant,1.1

    Parameters
    ----------
    groups_file_path: str
        output path for groups file
    wt_basenames
        the baseline ids to write to file
    mut_basenames
        the mutant ids to write to file

    """

    def process_specimens(specimen_list, genotype, crls):
        for volname in specimen_list:

            if use_crl:
                try:
                    extension_stipped_fname = common.strip_img_extension(volname)
                    staging_metric = crls.loc[extension_stipped_fname].value

                except KeyError as e:
                    logging.exception(f"Cannot find {volname} in the staging info file. {e}")
                    raise ValueError("Cannot find {} in the staging info file".format(volname))

                entries.append([volname, genotype, staging_metric])

            else:
                entries.append([volname, genotype])

    if wt_staging is not None and mut_staging is not None:
        use_crl = True
        columns = ['volume_id' ,'genotype', 'crl']


    else:
        use_crl = False
        columns = ['volume_id', 'genotype']

    entries = []

    # drop an WTs from mutant line if present in baseline set
    for bn in mut_basenames:
        if bn in wt_basenames:
            wt_basenames.remove(bn)

    process_specimens(wt_basenames, 'wildtype', wt_staging)
    process_specimens(mut_basenames, 'mutant', mut_staging)


    df = pd.DataFrame.from_records(entries, columns=columns)

    try:
        df.to_csv(groups_file_path, index=False)

    except FileNotFoundError as e:
        logging.exception("Cannot open combined groups file:\n".format(groups_file_path, e.strerror))
        raise LamaDataException("Cannot open combined groups file:\n".format(groups_file_path, e.strerror))


def get_formulas(config):
    """
    Extract the linear model formula from the stats config file. Just extract the independent varibale names for now

    example config entry = formulas: ['data ~ genotype']

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
    else:
        for formula_string in formulas:
            formula_elements = formula_string.split()[0::2][
                               1:]  # extract all the effects, miss out the dependent variable
            parsed_formulas.append(','.join(formula_elements))
        return parsed_formulas


def get_normalisation(config, mask_array):
    normalisation_roi = config.get('normalisation', 'mask')
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
                raise LamaDataException("Cannot read roi mask image for normalisation {}".format(n))
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
        raise IOError("Can't find mask {}. A mask is needed for the stats analysis".format(fixed_mask))

    # Look for inverted masks and add path to config if present
    # inverted_mask_dir

    voxel_size = config.get('voxel_size')
    if not voxel_size:
        voxel_size = common.DEFAULT_VOXEL_SIZE
        logging.warning("Voxel size not set in config. Using a default of {}".format(common.DEFAULT_VOXEL_SIZE))

    config.voxel_size = float(voxel_size)

    config.project_name = config.get('project_name', '_')

    config.mask_array_flat = common.img_path_to_array(fixed_mask).ravel()

    invert_config = config.get('invert_config_file')
    if invert_config:
        config.invert_config_file = join(root_dir, invert_config)

    # invert_mask = config.get('invert_config_file')
    # if invert_config:
    #     config.invert_config_file = join(root_dir, invert_config)

    config.blur_fwhm = config.get('blur_fwhm', DEFULAT_BLUR_FWHM)

    return config


def get_labels_and_names(root_dir: Path, label_map_path: Path, label_names_path: Path):
    """
    
    Parameters
    ----------
    root_dir: str
        
    label_map_path
    label_names_path

    Returns
    -------
    
    Notes
    -----
    see common.load_label_map_names for details on input file format

    """
    # Get the label maps and organ names, if used

    if label_map_path:
        lp = (root_dir / label_map_path).resolve()
        label_map = common.img_path_to_array(lp)
    else:
        label_map = None

    if label_names_path:
        label_names_path = (root_dir / label_names_path).resolve()
        organ_names = common.load_label_map_names(label_names_path)
    else:
        organ_names = None

    return label_map, organ_names


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

    stats_object = stats_method(outdir, analysis_prefix, config, analysis_config)

    # Run each dataset found in the stats.yaml config and rtun it against the appropraiate test
    for test in stats_tests:
        if test == 'LM' and not common.is_r_installed():
            logging.warning("Could not do linear model test for {}. Do you need to install R?".format(analysis_name))
            continue
        stats_object.run(STATS_METHODS[test], analysis_name)
        if config.invert_config_path:
            stats_object.invert(config.invert_config_path)
    del stats_object


