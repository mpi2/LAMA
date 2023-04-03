"""
The entry point to running permutation-based statistics.


Usage
-----

run_permutation_stats.run(test_data/registration_test_data/baseline,
                        data/registration_test_data/mutant,
                        test_data/stats_test_data/test_output,
                        1000)


Currently, this module works only with organ volume data. The voxel-based methods are currently too big to do this.
Think about including voxel-based stats in the future

Outline of pipeline
-------------------
Before running the permutation statistics we need to have run jobrunner.py on the baseline and mutant data.

The main function in this module run() calles the following functions during the pipeline:

get_organ_volume_data and get_staging_data
    search the registration output folders for the CSVs that contain the organ volumes
    and staging data and collate into single csvs.


distributions.null and distributions.alternative
    Use the dataframes from the precedding functions to generate null and alternative p-value distributiuon dataframes

p_thresholds.get_thresholds
    Using the null and alternative distributions, these functions generate organ-spceific p-value thresholds.
    These are generated for both line-level and specimen level calls.

annotate
    This function creates final results CSVs.
        Puts a line-level csv in the line/output/stats_/
        Puts specimen-level csv files in line/output/stats_/specimen_level

"""

from pathlib import Path
from datetime import date
import re
import pandas as pd
import numpy as np
from scipy.stats import zmap
from logzero import logger as logging
import yaml
from itertools import compress
from lama import common
from lama.stats.permutation_stats import distributions
from lama.stats.permutation_stats import p_thresholds
from lama.paths import specimen_iterator, get_specimen_dirs, LamaSpecimenData
from lama.qc.organ_vol_plots import make_plots, pvalue_dist_plots
from lama.common import write_array, read_array, init_logging, git_log, LamaDataException
from lama.stats.common import cohens_d
from lama.stats.penetrence_expressivity_plots import heatmaps_for_permutation_stats

GENOTYPE_P_COL_NAME = 'genotype_effect_p_value'
TREAT_P_COL_NAME = 'treatment_effect_p_value'
INTER_P_COL_NAME = 'interaction_effect_p_value'

PERM_SIGNIFICANT_COL_NAME = 'significant_cal_p'

PERM_SIGNIFICANT_COL_LIST = ['significant_cal_p_geno', 'significant_cal_p_treat', 'significant_cal_p_inter']

PERM_T_COL_NAME = 't'


def write_specimen_info(wt_wev, mut_wev, outfile):
    """
    Write a csv with some summary info on specimens
    currently only returns Z-score of mutants
    """

    def sortwev(x):
        return x
    wev_z = zmap(mut_wev.staging, wt_wev.staging)
    mut_wev['WEV_zscore'] = wev_z
    mut_wev.sort_values('WEV_zscore', key=sortwev, inplace=True)
    mut_wev.to_csv(outfile)


def get_radiomics_data(rad_dir: Path, wt_dir: Path, mut_dir: Path, treat_dir: Path, inter_dir: Path) -> pd.DataFrame:
    """
    Given a root registration directory, collate all the organ volume CSVs into one file.
    Write out the combined organ volume CSV into the root registration directory.

    Parameters
    ----------
    root_dir
        The path to the root registration directory

    Returns
    -------
    The combined data frame of all the organ volumes
    specimen id in index organs in rows
    """



    # get the features_per_embryo and convert it into per organs
    feature_dir = rad_dir / "features"
    #common.gather_rad_data(feature_dir)

    org_dir = rad_dir / "organs"

    # get the organ data and load it as one massive file
    file_names = [spec for spec in common.get_file_paths(folder=org_dir, extension_tuple=".csv") if "0." in str(spec)]
    file_names.sort()
    df_list = []


    staging = pd.concat([get_staging_data(_dir) for _dir in [wt_dir, mut_dir, treat_dir, inter_dir]]).rename(
        columns={'value': 'staging'})

    for org_name in file_names:
        # read dataset
        d = pd.read_csv(org_name, index_col=0).dropna(axis=1)
        # tag the columns with the organ_number

        # For some reason,  "." stuffs up the pipeline and adds a space, just remove it
        org = str(d.org[0]).replace(".0","")
        d.drop(columns=['HPE', 'genotype', 'background', 'org'], inplace=True)

        # patsy has a fit with "-" thinks I'm subtracting
        # I use '__' as a method to identifiy radiomics data
        d.set_axis([(col + '__' + org).replace("-","_") for col in d.columns], axis=1, inplace=True)
        d = d.reindex(staging.index)
        d.divide(staging['staging'])
        df_list.append(d)

    # horizontal merge - hope it works
    data = pd.concat(df_list, axis=1)
    #data = data.loc[:, data.columns.str.contains('shape')]
    data = data.loc[:, ~data.columns.str.contains('2D')]
    data = pd.concat([data, staging], axis=1)
    return data


def get_organ_volume_data(root_dir: Path) -> pd.DataFrame:
    """
    Given a root registration directory, collate all the organ volume CSVs into one file.
    Write out the combined organ volume CSV into the root registration directory.

    Parameters
    ----------
    root_dir
        The path to the root registration directory

    Returns
    -------
    The combined data frame of all the organ volumes
    specimen id in index organs in rows
    """

    dataframes = []

    s: LamaSpecimenData

    for s in get_specimen_dirs(root_dir):
        # for line_dir, specimen_dir in specimen_iterator(output_dir):

        # organ_vol_file = specimen_dir / 'output' / common.ORGAN_VOLUME_CSV_FILE
        organ_vol_file = s.outroot / common.ORGAN_VOLUME_CSV_FILE

        if not organ_vol_file.is_file():
            raise FileNotFoundError(f'Cannot find organ volume file {organ_vol_file}')

        df = pd.read_csv(organ_vol_file, index_col=0)

        if len(df) == 0:
            raise ValueError(f'{organ_vol_file} is empty')

        dataframes.append(df)

    # Write the concatenated organ vol file to single csv
    all_organs = pd.concat(dataframes)

    # outpath = output_dir / common.ORGAN_VOLUME_CSV_FILE
    # all_organs.to_csv(outpath)

    return all_organs


def get_staging_data(root_dir: Path) -> pd.DataFrame:
    """
    Given a root registration directory, collate all the staging CSVs into one file.
    Write out the combined organ volume CSV into the root registration directory.

    Parameters
    ----------
    root_dir
        The path to the root registration directory

    Returns
    -------
    The combined dataframe of all the organ volumes
    """
    # output_dir = root_dir / 'output'

    dataframes = []
    s: LamaSpecimenData

    for s in get_specimen_dirs(root_dir):
        # for line_dir, specimen_dir in specimen_iterator(output_dir):

        # staging_info = specimen_dir / 'output' / common.STAGING_INFO_FILENAME
        staging_info = s.outroot / common.STAGING_INFO_FILENAME

        if not staging_info.is_file():
            raise FileNotFoundError(f'Cannot find staging info file {staging_info}')

        df = pd.read_csv(staging_info, index_col=0)
        # df['line'] = line_dir.name
        df['line'] = s.line_id
        dataframes.append(df)

    # Write the concatenated staging info to the
    all_staging = pd.concat(dataframes)
    # outpath = output_dir / common.STAGING_INFO_FILENAME
    outpath = root_dir / common.STAGING_INFO_FILENAME
    all_staging.to_csv(outpath)

    return all_staging


def annotate(thresholds: pd.DataFrame,
             lm_results: pd.DataFrame,
             lines_root_dir: Path,
             is_line_level: bool = True,
             label_info: Path = None,
             label_map: Path = None,
             write_thresholded_inv_labels=False,
             fdr_threshold: float = 0.05,
             t_values: pd.DataFrame = None,
             organ_volumes: pd.DataFrame = None,
             two_way: bool = False,
             main_of_two_way: bool = False) -> pd.DataFrame:
    """
    Using the p_value thresholds and the linear model p-value results,
    create the following CSV files

        Line-level results
        specimen-level results

    Parameters
    ----------
    thresholds
        columns label(index), p_thresh, fdr, num_hits_across_all_lines/specimens
    lm_results
        The alternative distribution
        index: line/specimen id
        cols: labels (+ line_id for specimen_level)
    lines_root_dir
        The root directory to save the annotated CSV files. Each line to go in a subfolder
    is_line_level
        if not True, place results in specimen-level sub directory
    label_info
        CSV to map label number to name
    t_values
         same format as lm_results but with t-statistics
    organ_volumes
        All the organ volumes for baselines and mutants (as it was used in lm(), so probably normalised to whole embryo

    Returns
    -------
    Aggregated hit dataframe

    Notes
    -----
    TODO: Add file number prefixes so we don't overwrite mulyiple analyses done on the same day
    TODO: the organ_volumes folder name is hard-coded. What about if we add a new analysis type to the  permutation stats pipeline?
    """
    hit_dataframes = []

    if label_map:
        label_map = read_array(label_map)

    if two_way:
        thresholds = thresholds.pivot(columns='effect')
    # Iterate over each line or specimen (for line or specimen-level analysis)
    for id_, row in lm_results.iterrows():

        # Create a dataframe containing a p-value column. each row an organ
        df = row.to_frame()
        if not is_line_level:
            # specimen-level has an extra line column we need to remove
            df = df.T.drop(columns=['line']).T

        # Rename the line_specimen column to be more informative

        if (two_way and not main_of_two_way):
            df.drop(labels=['line'], axis=0, errors='ignore', inplace=True)

            try:
                fixed_vals = np.stack(df[id_])
                df['genotype_effect_p_value'] = pd.to_numeric(fixed_vals[:, 0])
                df['interaction_effect_p_value'] = pd.to_numeric(fixed_vals[:, 2])
                df['treatment_effect_p_value'] = pd.to_numeric(fixed_vals[:, 1])
                df.drop(columns=['two_way'], errors='ignore', inplace=True)
                df.drop(labels=['line'], axis=0, errors='ignore', inplace=True)

            except IndexError:
                # data wrangling - remove brackets and convert values to float
                fixed_vals = pd.DataFrame([re.sub('\[|\]', '', val).split(' ')[0:3]
                                           for val in df[id_]], index=df.index)
                df['genotype_effect_p_value'] = pd.to_numeric(fixed_vals[0], errors='coerce')
                df['interaction_effect_p_value'] = pd.to_numeric(fixed_vals[2], errors='coerce')
                df['treatment_effect_p_value'] = pd.to_numeric(fixed_vals[1], errors='coerce')

        # fix up the specimen main two-ways
        elif main_of_two_way:
            df.drop(labels=['line'], axis=0, errors='ignore', inplace=True)

            try:
                spec_name = df.columns
                df = pd.DataFrame(np.stack(df.iloc[:, 0]), index=df.index)
                # print("fixed_val ", fixed_vals, type(fixed_vals))
                # df = pd.DataFrame(fixed_vals, index=df.index)
                # print("numeric val", df)
                df.rename(columns={0: GENOTYPE_P_COL_NAME}, inplace=True)

            except IndexError:
                # this is only really for testing where the the arrays are not properly written by to_csv

                fixed_vals = pd.DataFrame([re.sub('\[|\]', '', val) for val in df.iloc[:, 0]], index=df.index)
                df = pd.DataFrame(pd.to_numeric(fixed_vals[0]), index=df.index)
                df.rename(columns={0: GENOTYPE_P_COL_NAME}, inplace=True)

        else:
            df.rename(columns={id_: GENOTYPE_P_COL_NAME}, inplace=True)

        if is_line_level:
            line = id_
        else:
            line = row['line']
            spec_id = id_

        # Merge the permutation results (p-thresh, fdr, number of hit lines for this label) with the mutant results

        df.index = df.index.astype(str)  # Index needs to be cast from object to enable merge
        df = df.merge(thresholds, left_index=True, right_index=True, validate='1:1')

        df.index.name = 'label'

        # Merge the t-statistics
        if t_values is not None:
            t_df = pd.DataFrame(t_values.loc[id_])

            if t_df.shape[1] != 1:  # We get multiple columns f there are duplicate specimen/line ids
                raise ValueError("Duplicate specimen names not allowed")

            t_df.columns = ['t']
            t_df.drop(columns=['line'], errors='ignore', inplace=True)  # this is for speciem-level results

            t_df.index = t_df.index.astype(str) # index must be string for radiomics data, organ data doesn't seem to care?

            df = df.merge(t_df, left_index=True, right_index=True, validate='1:1')
            if len(df) < 1:
                logging.info(f'skipping {id_} no hits')  # Should we continue at this point?

        # Add mean organ vol difference and cohens d
        df['mean_vol_ratio'] = None
        if is_line_level:
            df['cohens_d'] = None

        for label, row in df.iterrows():

            # Organ vols are prefixed with x so it can work with statsmodels
            label_col = f'{label}'if str(label).__contains__("__") else f'x{label}'
            label_organ_vol = organ_volumes[[label_col, 'line']]

            wt_ovs = label_organ_vol.loc[label_organ_vol.line == 'baseline',label_col]

            if two_way or main_of_two_way:
                # I think this is the only way to get the combs....

                # Giving ChatGPT a chance to shine - it loves using .loc, checking data is not null convert it to numpy
                mut_ovs = label_organ_vol.loc[
                    (label_organ_vol.line == 'mutants') & label_organ_vol[label_col].notnull(),
                    label_col
                ].to_numpy()

                treat_ovs = label_organ_vol.loc[
                    (label_organ_vol.line == 'treatment') & label_organ_vol[label_col].notnull(),
                    label_col
                ].to_numpy()

                int_ovs = label_organ_vol.loc[
                    (label_organ_vol.line == 'mut_treat') & label_organ_vol[label_col].notnull(),
                    label_col
                ].to_numpy()

                non_int_ovs = label_organ_vol.loc[
                    label_organ_vol.line.isin(['baseline', 'mutants', 'treatment']) & label_organ_vol[
                        label_col].notnull(),
                    label_col
                ].to_numpy()

                #line.values should be string
                if 'mut_treat' in label_organ_vol.line.values:
                    num_ovs = int_ovs
                    dem_ovs = non_int_ovs
                elif 'treatment' in label_organ_vol.line.values:
                    num_ovs = treat_ovs
                    dem_ovs = wt_ovs
                else:
                    num_ovs = mut_ovs
                    dem_ovs = wt_ovs

                # Specimen level - overwrite the num_ovs to be the single emb of interest
                if not is_line_level and two_way:
                    num_ovs = label_organ_vol.loc[label_organ_vol.index == row.index[0], label_col]

                elif not is_line_level and main_of_two_way:
                    num_ovs = label_organ_vol.loc[label_organ_vol.index == spec_name[0], label_col]


                df.loc[label, 'mean_vol_ratio'] = num_ovs.mean() / dem_ovs.mean()
                if is_line_level:
                    df.loc[label, 'cohens_d'] = cohens_d(num_ovs, dem_ovs)

            else:
                mut_ovs = label_organ_vol[label_organ_vol.line == line][label_col]

                df.loc[label, 'mean_vol_ratio'] = mut_ovs.mean() / wt_ovs.mean()
                if is_line_level:
                    cd = cohens_d(mut_ovs, wt_ovs)
                    df.loc[label, 'cohens_d'] = cd

        output_name = f'{id_}_organ_volumes_{str(date.today())}.csv'

        line_output_dir = lines_root_dir / line
        line_output_dir.mkdir(exist_ok=True)

        if not is_line_level:
            # If dealing with specimen-level stats, make subfolder to put results in
            line_output_dir = line_output_dir / 'specimen_level' / id_
            line_output_dir.mkdir(parents=True, exist_ok=True)

        output_path = line_output_dir / output_name


        add_two_way_significance(df, fdr_threshold) if two_way else add_significance(df, fdr_threshold)

        if label_info:
            df = add_label_names(df, label_info)

        df.to_csv(output_path)

        if two_way:
            # print(any(df[PERM_SIGNIFICANT_COL_LIST] == True, axis = 'columns'))

            eff_there = [(GENOTYPE_P_COL_NAME in df.columns),
                         (TREAT_P_COL_NAME in df.columns),
                         (INTER_P_COL_NAME in df.columns)]

            PERM_COL_LIST = list(compress(PERM_SIGNIFICANT_COL_LIST, eff_there))

            hit_df = df[(df[PERM_COL_LIST] == True).any(axis='columns')]
        else:
            hit_df = df[df['significant_cal_p'] == True]
            hit_df['line'] = line

        if not is_line_level:
            hit_df['specimen'] = spec_id

        hit_dataframes.append(hit_df)

        hit_labels_out = line_output_dir / f'{line}__hit_labels.nrrd'

        hits = hit_df.index

        if write_thresholded_inv_labels:
            _write_thresholded_label_map(label_map, hits, hit_labels_out)

    collated_df = pd.concat(hit_dataframes)
    return collated_df


def _write_thresholded_label_map(label_map: np.ndarray, hits, out: Path):
    """
    Write a label map with only the 'hit' organs in it
    """
    if label_map is None:
        return

    if len(hits) > 0:
        # Make a copy as it may be being used elsewhere
        l = np.copy(label_map)
        # Clear any non-hits
        l[~np.isin(l, hits)] = 0

        write_array(l, out)


def add_label_names(df: pd.DataFrame, label_info: Path) -> pd.DataFrame:
    """
    Added label names to hits dataframe with merge on label metadata
    """
    label_df = pd.read_csv(label_info, index_col=0)
    #if its radiomics data, the columns will have __
    if df.index[0].__contains__("__"):# this is for radiomics data
        # 3D stuffs up labelling
        label_nums = [int(re.findall('\d+', _row.replace('3D', ""))[0]) for _row in df.index]
        df['label_name'] = [label_df.loc[num]['label_name'] for num in label_nums]
        # so this just adds the label_name and no_analysis columns by matching the label number with the feature
        if 'no_analysis' in label_df:
            df['no_analysis'] = [label_df.loc[num]['no_analysis'] for num in label_nums]
    else:
        label_df.index = label_df.index.astype(str)
        df = df.merge(right=label_df[['label_name']], left_index=True, right_index=True)
        if 'no_analysis' in label_df:
            df = df.merge(right=label_df[['no_analysis']], left_index=True, right_index=True)
    return df


def add_significance(df: pd.DataFrame, threshold: float):
    """
    Add a significance column to the output csv in place.
    Set significance to True if the genotype p-value is lower than the p threshold for that organ
    and the fdr is lower than fdr threshold.
    And sort values by significance
    """
    df[PERM_SIGNIFICANT_COL_NAME] = (df[GENOTYPE_P_COL_NAME] <= df['p_thresh']) & (df['fdr'] <= threshold)

    df.sort_values(by=[PERM_SIGNIFICANT_COL_NAME, GENOTYPE_P_COL_NAME], ascending=[False, True], inplace=True)


def add_two_way_significance(df: pd.DataFrame, threshold: float):
    """
    Add a significance column to the output csv in place.
    Set significance to True if the genotype p-value is lower than the p threshold for that organ
    and the fdr is lower than fdr threshold.
    And sort values by significance
    """

    eff_there = [(GENOTYPE_P_COL_NAME in df.columns),
                 (TREAT_P_COL_NAME in df.columns),
                 (INTER_P_COL_NAME in df.columns)]

    P_COL_LIST = [('genotype', GENOTYPE_P_COL_NAME),
                  ('treatment', TREAT_P_COL_NAME),
                  ('interaction', INTER_P_COL_NAME)]

    # cond_list = ['genotype', 'treatment', 'interaction']
    sort_list = list(compress([False, False, False], eff_there))

    PERM_COL_LIST = list(compress(PERM_SIGNIFICANT_COL_LIST, eff_there))

    for i, cond in enumerate(list(compress(P_COL_LIST, eff_there))):
        df[PERM_COL_LIST[i]] = (df[cond[1]] <= df[('p_thresh', cond[0])]) \
                               & (df[('fdr', cond[0])] <= threshold)

    df.sort_values(by=PERM_COL_LIST, ascending=sort_list, inplace=True)


def prepare_data(wt_organ_vol: pd.DataFrame,
                 wt_staging: pd.DataFrame,
                 mut_organ_vol: pd.DataFrame,
                 mut_staging: pd.DataFrame,
                 label_meta: Path = None,
                 normalise_to_whole_embryo=False,
                 qc_file: Path = None,
                 two_way_data: list = [],
                 two_way: bool = False,
                 rad_data: bool = False) -> pd.DataFrame:
    """
    Merge the mutant and wildtype dtaframes
    Optionally normalise to staging metric (Usually whole embryo volume)
    Optionally remove any qc-flagged organs (These will be set to 'nan')

    Returns
    -------
    Dataframe with following columns:
        - a column for each label (prefixed with 'x' as statsmodels does not like integer ids)
        - line
        - genotype (baseline or mutant)
        - staging (whole embryo volume)
    """

    wt_staging.rename(columns={'value': 'staging'}, inplace=True)
    mut_staging.rename(columns={'value': 'staging'}, inplace=True)
    wt_staging.index = wt_staging.index.astype(str)

    if two_way:
        # unpack data
        treat_staging, treat_organ_vol, inter_staging, inter_organ_vol = two_way_data
        # Now do essentially the same stuff as wt and muts
        treat_staging.rename(columns={'value': 'staging'}, inplace=True)
        inter_staging.rename(columns={'value': 'staging'}, inplace=True)

    else:
        # just the one-way
        # Ensure all indices are same type
        for d in [wt_organ_vol, mut_organ_vol, wt_staging, mut_staging]:
            d.index = d.index.astype(str)

    if normalise_to_whole_embryo:
        logging.info('Normalising organ volume to whole embryo volume')
        wt_organ_vol = wt_organ_vol.divide(wt_staging['staging'], axis=0)
        mut_organ_vol = mut_organ_vol.divide(mut_staging['staging'], axis=0)
        if two_way:
            # normalise the other groups
            treat_organ_vol = treat_organ_vol.divide(treat_staging['staging'], axis=0)
            inter_organ_vol = inter_organ_vol.divide(inter_staging['staging'], axis=0)

    # merge the organ vol
    # list comprehension cause why not
    organ_vols = pd.concat([wt_organ_vol, mut_organ_vol, treat_organ_vol, inter_organ_vol]) if two_way \
        else pd.concat([wt_organ_vol, mut_organ_vol])

    # Drop any organ columns that has only zero values. These are the gaps in the label map caused by merging labels
    # in the atlas
    organ_vols = organ_vols.loc[:, (organ_vols != 0).any(axis=0)]

    # For the statsmodels linear mode to work, column names cannot start with a digit. Prefix with 'x'
    organ_vols.columns = [f'x{x}' if x.isdigit() else x for x in organ_vols.columns]

    staging = pd.concat([wt_staging, mut_staging, treat_staging, inter_staging]) if two_way \
        else pd.concat([wt_staging, mut_staging])

    # Merge staging to the organvolume dataframe. First drop line so we don't get duplicate entries
    # staging.drop(columns=['line'], inplace=True)

    data = pd.concat([organ_vols, staging], axis=1)

    # Filter any labels that have been flagged at the label-level (for all specimens)
    if label_meta:

        label_meta = pd.read_csv(label_meta, index_col=0)

        if 'no_analysis' in label_meta:
            # If we have a no_analysis column, drop labels that are flagged

            flagged_labels = label_meta[label_meta.no_analysis == True].index

            data.drop(columns=[f'x{x}' for x in flagged_labels if f'x{x}' in data], inplace=True)

    # QC-flagged organs from specimens specified in QC file are set to None
    if qc_file:
        logging.info(f'Excluding organ volumes specified in: {qc_file}')
        qc = pd.read_csv(qc_file)

        for _, row in qc.iterrows():
            qc_id = str(row.id)

            if qc_id not in data.index:
                raise LamaDataException(f'QC flagged specimen {row.id} does not exist in dataset')

            if f'x{row.label}' not in data:
                raise LamaDataException(f'QC flagged label, {row.label}, does not exist in dataset')

            data.loc[qc_id, f'x{row.label}'] = None

    return data


def run(wt_dir: Path,
        mut_dir: Path,
        out_dir: Path,
        num_perms: int,
        label_info: Path = None,
        label_map_path: Path = None,
        line_fdr: float = 0.05,
        specimen_fdr: float = 0.2,
        normalise_to_whole_embryo: bool = True,
        qc_file: Path = None,
        voxel_size: float = 1.0,
        two_way: bool = False,
        treat_dir: Path = None,
        inter_dir: Path = None,
        rad_dir: Path = None):
    """
    Run the permutation-based stats pipeline

    Parameters
    ----------
    wt_dir
        Root of the wild type registration output
        This should contain an 'output' folder that contains a single baseline folder that contains multiple specimen folders
    mut_dir
        Root of the mutant registration output
        This should contain 'output' folder that contains multiple mutant lines folder, each containing one or more mutant specimen folders
    out_dir
        Where to store the intermediate results of the permutation testing
    num_perms
        number of permutations to do
    log_dependent
        if True, apply numpy.log to all the dependent values (organ volumes)
    label_info
        if supplied, use it to annotate the results with label names. Also can be used to filter certain labels from the
        analysis using the 'no_analysis' column
    line_fdr
        the FDR threshold at which to accept line level calls
    specimen_fdr
        the FDR threshold at which to accept specimen-level calls
    normalise_to_whole_embryo:
        Whether to divide the organ each organ volume by whole embryo volume
    qc_file
        csv indicating labels from specimens that should be excluded from the analysis
        columns:
        - id: the specimen id
        - line: the line id
        - label: the label to exclude (int)
        - label_name (optional)
    voxel_size
        For calcualting organ volumes
    two_way
        Activates the two-way simulation
    """
    # Collate all the staging and organ volume data into csvs
    np.random.seed(999)
    init_logging(out_dir / 'stats.log')
    logging.info(git_log())
    logging.info(f'Running {__name__} with following commands\n{common.command_line_agrs()}')


    # data
    # index: spec_id
    # cols: label_nums, with staging and line columns at the end
    if rad_dir:
        logging.info('Searching for staging data')
        wt_staging = get_staging_data(wt_dir)
        mut_staging = get_staging_data(mut_dir)
        logging.info('Collecting Radiomics data')
        data = get_radiomics_data(rad_dir, wt_dir, mut_dir, treat_dir, inter_dir)
        # turn on textures at your own risk
        data.to_csv(out_dir / 'radiomics_data.csv')

    else:
        logging.info('Searching for staging data')
        wt_staging = get_staging_data(wt_dir)
        mut_staging = get_staging_data(mut_dir)

        logging.info('searching for organ volume data')
        wt_organ_vol = get_organ_volume_data(wt_dir)
        mut_organ_vol = get_organ_volume_data(mut_dir)
        if two_way:
            logging.info('Searching for two-way staging and organ volume data')
            treat_staging = get_staging_data(treat_dir)
            inter_staging = get_staging_data(inter_dir)
            treat_organ_vol = get_organ_volume_data(treat_dir)
            inter_organ_vol = get_organ_volume_data(inter_dir)
            two_way_data = [treat_staging, treat_organ_vol,
                            inter_staging, inter_organ_vol]
        else:
            two_way_data = []
        data = prepare_data(wt_organ_vol,
                        wt_staging,
                        mut_organ_vol,
                        mut_staging,
                        label_meta=label_info,
                        normalise_to_whole_embryo=normalise_to_whole_embryo,
                        qc_file=qc_file,
                        two_way=two_way,
                        two_way_data=two_way_data)

        data.to_csv(out_dir / 'input_data.csv')



    # get rad data



    # Make plots
    # data_for_plots = data.copy()
    # data_for_plots.columns = [x.strip('x') for x in data_for_plots.columns]  # Strip any xs
    # # If data has been normalised to WEV revert back for plots
    # if normalise_to_whole_embryo:
    #     for col in data_for_plots.columns:
    #         if col.isdigit():
    #             data_for_plots[col] = data_for_plots[col] * data_for_plots['staging']
    lines_root_dir = out_dir / 'lines'
    # make_plots(data_for_plots, label_info, lines_root_dir, voxel_size=voxel_size)

    # Keep a record of the input data used in the analsysis


    # Keep raw data for plotting
    # raw_wt_vols = wt_organ_vol.copy()   # These includes QCd speciemns need to remove

    out_dir.mkdir(exist_ok=True, parents=True)  # Root directory for output

    # make directory to store distributions and thresholds
    dists_out = out_dir / 'distributions'
    dists_out.mkdir(exist_ok=True)

    # Get the null distributions
    logging.info('Generating null distribution')
    line_null, specimen_null = distributions.null(data, num_perms, two_way=two_way)

    # with open(dists_out / 'null_ids.yaml', 'w') as fh:
    #     yaml.dump(null_ids, fh)

    null_specimen_pvals_file = dists_out / 'null_specimen_dist_pvalues.csv'

    null_line_pvals_file = dists_out / 'null_line_dist_pvalues.csv'

    # Write the null distributions to file
    line_null.to_csv(null_line_pvals_file)
    specimen_null.to_csv(null_specimen_pvals_file)

    # Get the alternative p-value distribution (and t-values now (2 and 3)
    logging.info('Generating alternative distribution')
    line_alt, spec_alt, line_alt_t, spec_alt_t = distributions.alternative(data, two_way=two_way)

    line_alt_pvals_file = dists_out / 'alt_line_dist_pvalues.csv'
    spec_alt_pvals_file = dists_out / 'alt_specimen_dist_pvalues.csv'

    # Write the alternative distributions to file
    line_alt.to_csv(line_alt_pvals_file)
    spec_alt.to_csv(spec_alt_pvals_file)

    line_organ_thresholds = p_thresholds.get_thresholds(line_null, line_alt, two_way=two_way)

    line_thresholds_path = dists_out / 'line_organ_p_thresholds.csv'
    line_organ_thresholds.to_csv(line_thresholds_path)

    # let's tidy up our data from the specimen calls in the two_way
    if two_way:
        # TODO: Don't hard-code this
        specimen_inter_nulls = specimen_null[specimen_null.iloc[:, 0].str.len() == 3]

        specimen_main_nulls = specimen_null[specimen_null.iloc[:, 0].str.len() == 1]
        specimen_geno_nulls, specimen_treat_nulls = np.vsplit(specimen_main_nulls, 2)


        specimen_inter_alt = spec_alt[spec_alt.iloc[:, 1].str.len() == 3]


        specimen_main_alt = spec_alt[spec_alt.iloc[:, 1].str.len() == 1]



        # so firstly let's get the names and conditions from the data
        group_info = data['line']


        # TODO: think whether to truly put mut_treat in main comparisons
        mut_names = group_info[(group_info == 'mutants') | (group_info == 'mut_treat')].index
        treat_names = group_info[(group_info == 'treatment') | (group_info == 'mut_treat')].index

        specimen_geno_alt = specimen_main_alt[specimen_main_alt.index.isin(mut_names)]
        specimen_treat_alt = specimen_main_alt[specimen_main_alt.index.isin(treat_names)]

        geno_alt_path = dists_out / 'specimen_geno_pvals.csv'
        treat_alt_path = dists_out / 'specimen_treat_pvals.csv'
        inter_alt_path = dists_out / 'specimen_inter_pvals.csv'

        specimen_geno_alt.to_csv(geno_alt_path)
        specimen_treat_alt.to_csv(treat_alt_path)
        specimen_inter_alt.to_csv(inter_alt_path)

        geno_thresholds = p_thresholds.get_thresholds(specimen_geno_nulls, specimen_geno_alt, two_way=two_way)
        treat_thresholds = p_thresholds.get_thresholds(specimen_treat_nulls, specimen_treat_alt, two_way=two_way)
        inter_thresholds = p_thresholds.get_thresholds(specimen_inter_nulls, specimen_inter_alt, two_way=two_way)

        geno_thresholds_path = dists_out / 'specimen_geno_p_thresholds.csv'
        treat_thresholds_path = dists_out / 'specimen_treat_p_thresholds.csv'
        inter_thresholds_path = dists_out / 'specimen_inter_p_thresholds.csv'

        geno_thresholds.to_csv(geno_thresholds_path)
        treat_thresholds.to_csv(treat_thresholds_path)
        inter_thresholds.to_csv(inter_thresholds_path)

    else:
        specimen_organ_thresholds = p_thresholds.get_thresholds(specimen_null, spec_alt, two_way=two_way)
        spec_thresholds_path = dists_out / 'specimen_organ_p_thresholds.csv'
        specimen_organ_thresholds.to_csv(spec_thresholds_path)

    logging.info('Annotating lines')

    # The lines root doesn't really exist in a two-way
    if two_way:
        lines_root_dir = out_dir
    else:
        lines_root_dir = out_dir / 'lines'
        lines_root_dir.mkdir(exist_ok=True)

    # Annotate lines
    logging.info(f"Annotating lines, using a FDR threshold of {line_fdr}")
    line_hits = annotate(line_organ_thresholds, line_alt, lines_root_dir, label_info=label_info,
                         label_map=label_map_path, write_thresholded_inv_labels=True, fdr_threshold=line_fdr,
                         t_values=line_alt_t,
                         organ_volumes=data, two_way=two_way)

    line_hits.to_csv(out_dir / 'line_hits.csv')

    # Annotate specimens
    logging.info(f"Annotating specimens, using a FDR threshold of {specimen_fdr}")
    if two_way:
        geno_spec_hits = annotate(geno_thresholds, specimen_geno_alt, lines_root_dir, is_line_level=False,
                                  label_info=label_info, label_map=label_map_path, fdr_threshold=specimen_fdr,
                                  t_values=spec_alt_t,
                                  organ_volumes=data, main_of_two_way=True)

        treat_spec_hits = annotate(treat_thresholds, specimen_treat_alt, lines_root_dir, is_line_level=False,
                                   label_info=label_info, label_map=label_map_path, fdr_threshold=specimen_fdr,
                                   t_values=spec_alt_t,
                                   organ_volumes=data, main_of_two_way=True)

        inter_spec_hits = annotate(inter_thresholds, specimen_inter_alt, lines_root_dir, is_line_level=False,
                                   label_info=label_info, label_map=label_map_path, fdr_threshold=specimen_fdr,
                                   t_values=spec_alt_t,
                                   organ_volumes=data, two_way=True)
        geno_spec_hits.to_csv(out_dir / 'specimen_level_geno_hits.csv')
        treat_spec_hits.to_csv(out_dir / 'specimen_level_treat_hits.csv')
        inter_spec_hits.to_csv(out_dir / 'specimen_level_inter_hits.csv')

    else:
        spec_hits = annotate(specimen_organ_thresholds, spec_alt, lines_root_dir, is_line_level=False,
                             label_info=label_info, label_map=label_map_path, fdr_threshold=specimen_fdr,
                             t_values=spec_alt_t,
                             organ_volumes=data)

        spec_hits.to_csv(out_dir / 'specimen_level_hits.csv')

    # Make plots
    data_for_plots = data.copy()
    data_for_plots.columns = [x.strip('x') for x in data_for_plots.columns]  # Strip any xs
    # If data has been normalised to WEV revert back for plots
    if normalise_to_whole_embryo:
        for col in data_for_plots.columns:
            if col.isdigit():
                data_for_plots[col] = data_for_plots[col] * data_for_plots['staging']

    make_plots(data_for_plots, label_info, lines_root_dir, voxel_size=voxel_size, two_way=two_way)

    # Get specimen info. Currently just the WEV z-score to highlight specimens that are too small/large
    spec_info_file = out_dir / 'specimen_info.csv'
    #write_specimen_info(wt_staging, mut_staging, spec_info_file)

    dist_plot_root = out_dir / 'distribution_plots'
    line_plot_dir = dist_plot_root / 'line_level'
    line_plot_dir.mkdir(parents=True, exist_ok=True)
    pvalue_dist_plots(line_null, line_alt, line_organ_thresholds, line_plot_dir, two_way=two_way)

    specimen_plot_dir = dist_plot_root / 'specimen_level'
    specimen_plot_dir.mkdir(parents=True, exist_ok=True)
    if two_way:
        # fix up vals.
        pvalue_dist_plots(specimen_geno_nulls, specimen_geno_alt.drop(columns=['line']), geno_thresholds,
                          specimen_plot_dir, main_of_two_way=True)
        pvalue_dist_plots(specimen_treat_nulls, specimen_treat_alt.drop(columns=['line']), treat_thresholds,
                          specimen_plot_dir, main_of_two_way=True)
        pvalue_dist_plots(specimen_inter_nulls, specimen_inter_alt.drop(columns=['line']), inter_thresholds,
                          specimen_plot_dir, two_way=True)
    else:
        pvalue_dist_plots(specimen_null, spec_alt.drop(columns=['line']), specimen_organ_thresholds, specimen_plot_dir)


    rad_plot = True if rad_dir else False
    heatmaps_for_permutation_stats(lines_root_dir, two_way=two_way, label_info_file=label_info, rad_plot=rad_plot)
