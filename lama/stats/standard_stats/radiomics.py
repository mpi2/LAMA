import nrrd
from lama import common
import os
import numpy as np
from radiomics import featureextractor
import SimpleITK as sitk
import pandas as pd
import logging
import scipy


def extract_registrations(root_dir, labs_of_interest=None):
    '''

    either extracts the rigid registrations (i.e. the volumes)
    or specific labels (if the labels are specified)

    Parameters
    ----------
    root_dir
    labels - label of interest

    Returns
    -------
    list of either sitk images (rigid regs) or direct labels as arrays (labels)

    '''
    if labs_of_interest:
        # extract the inverted labels of interest
        label_paths = [spec_path for spec_path in common.get_file_paths(root_dir) if
                       ('inverted_labels' in str(spec_path))]

        label_paths.sort(key=lambda x: os.path.basename(x))

        # empty list
        extracts = [None] * len(label_paths)

        for i, path in enumerate(label_paths):
            # clean label_files to only contain orgs of interest
            label = common.LoadImage(path).img
            label_arr = sitk.GetArrayFromImage(label)
            # I think its better to just grab the single files for all orgs
            # then separate the labels during radiomics calculations
            label_arr[~np.isin(label_arr, labs_of_interest)] = 0
            extracts[i] = label_arr

    else:
        # extract the rigid
        reg_paths = [spec_path for spec_path in common.get_file_paths(root_dir) if ('registrations' in str(spec_path))]
        rigid_paths = [spec_path for spec_path in reg_paths if ('rigid' in str(spec_path))]
        rigid_paths.sort(key=lambda x: os.path.basename(x))

        # just an easy way to load the images
        extracts = [common.LoadImage(path) for path in rigid_paths]

    return extracts


def calc_all_features(target_dir, labs_of_int=None):
    '''i
    Performs the pyradiomic calculations


    Parameters
    ----------
    target_dir

    labs_of_int

    Returns
    -------

    '''

    # fix label input

    if labs_of_int != None:
        labs_of_int = [float(i) for i in labs_of_int.split(",")] if "," in labs_of_int else [float(labs_of_int)]
    else:
        labs_of_int = list(range(1, 210))

        # extract rigid registrations and inverted labels
    logging.info("Extracting Rigids")
    rigids = extract_registrations(target_dir)

    logging.info("Extracting Inverted Labels")
    labels = extract_registrations(target_dir, labs_of_int)

    # Get the radiomic measurements


    full_results = pd.Series([])
    names = [x.img_path for x in rigids]
    images = [x.img for x in rigids]

    # TODO: reduce dimensionality?
    for i, org in enumerate(labs_of_int):
        org_results = pd.Series([])
        for j, spec_label in enumerate(labels):
            # remove other labels
            arr_spec = np.copy(spec_label)
            arr_spec[spec_label != org] = 0
            arr_spec[spec_label == org] = 1

            if np.count_nonzero(arr_spec) < 1000:
                print("null label")
                continue


            mask = sitk.GetImageFromArray(arr_spec)

            # make sure its in the same orientation as the image
            mask.CopyInformation(images[j])



            extractor = featureextractor.RadiomicsFeatureExtractor()
            extractor.enableAllImageTypes()
            extractor.enableAllFeatures()
            result = extractor.execute(images[j], mask)

            features = pd.DataFrame.from_dict(result, orient='index',
                                              columns=[str(os.path.splitext(os.path.basename(names[j]))[0])])

            # transpose so features are columns
            features = features.transpose()

            # remove diagnostic columns and add
            features = features[features.columns.drop(list(features.filter(regex="diagnostics")))].add_suffix(('_'+str(org)))

            #org_results.append(features)
            org_results = pd.concat([org_results, features], axis=0)


        full_results=pd.concat([full_results, org_results], axis=1)
        print(full_results)

    # _metadata = features.index.str.split('_', expand=True).to_frame(index=False,
    #                                                                name=['Date', 'Exp', 'Contour_Method',
    #                                                                      'Tumour_Model', 'Position', 'Age',
    #                                                                      'Cage_No.', 'Animal_No.'])
    # _metadata.reset_index(inplace=True, drop=True)
    # features.reset_index(inplace=True, drop=True)
    # features = pd.concat([_metadata, features], axis=1)

    # features.index.rename('scanID', inplace=True)

    return features



