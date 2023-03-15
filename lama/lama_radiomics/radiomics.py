import nrrd
from lama import common
import os
import numpy as np
from radiomics import featureextractor
import SimpleITK as sitk
import pandas as pd
from logzero import logger as logging
from pathlib import Path
from filelock import SoftFileLock, Timeout
import socket
from datetime import datetime
import sys
import signal
import tempfile
from lama.monitor_memory import MonitorMemory
from lama.img_processing import normalise
from scipy import ndimage
import raster_geometry as rg
import subprocess


JOBFILE_NAME = 'radiomics_jobs.csv'


def extract_registrations(root_dir, labs_of_interest=None, norm_label=None,  fnames = None, stats_mask: bool=False):
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
    rad_dir = root_dir / "radiomics_output"
    os.makedirs(rad_dir, exist_ok=True)

    if labs_of_interest:
        # save to label folder

        outdir_string = "stage_labels" if norm_label else "stats_mask" if stats_mask else "inverted_labels"
        query_string = 'inverted_stats_mask' if stats_mask else outdir_string

        outdir = rad_dir / outdir_string
        os.makedirs(outdir, exist_ok=True)

        # extract the inverted labels of interest
        file_paths = [spec_path for spec_path in common.get_file_paths(root_dir) if
                      (query_string in str(spec_path))]

        logging.info(rad_dir)

        #file_paths.sort(key=lambda x: os.path.basename(x))

        # empty list
        with tempfile.NamedTemporaryFile() as ntf:
            temp_name = ntf.name
            extracts = np.memmap(ntf,shape=(len(file_paths),), dtype=object)


        for i, path in enumerate(file_paths):
            # clean label_files to only contain orgs of interest
            label = common.LoadImage(path).img

            label_arr = sitk.GetArrayFromImage(label)
            t = tempfile.TemporaryFile()
            m = np.memmap(t, dtype=label_arr.dtype, mode='w+', shape=label_arr.shape)
            m[:] = label_arr

            # I think its better to just grab the single files for all orgs
            # then separate the labels during radiomics calculations
            if not stats_mask:
                m[~np.isin(label_arr, labs_of_interest)] = 0
            extracts[i] = sitk.GetImageFromArray(m)
            extracts[i].CopyInformation(label)

    else:
        # extract the rigid
        outdir = rad_dir / "rigids"
        os.mkdir(outdir)

        reg_paths = [spec_path for spec_path in common.get_file_paths(root_dir) if ('registrations' in str(spec_path))]
        file_paths = [spec_path for spec_path in reg_paths if ('rigid' in str(spec_path))]


        # just an easy way to load the images
        extracts = [common.LoadImage(path).img for path in file_paths]

    #sort file paths
    #file_paths.sort(key=lambda x: os.path.basename(x))
    # write to new_folder for job file / increase debugging speed
    for i, vol in enumerate(extracts):

        file_name = str(Path(outdir / os.path.basename(file_paths[i])))

        #print("vol.img", vol.img)
        logging.info("Writing : {}". format(file_name))
        sitk.WriteImage(vol, file_name, useCompression=True)


    return extracts


def make_rad_jobs_file(jobs_file: Path, file_paths: list):
    """
    Creates a joblist csv file for use with the radiomics pipeline.
    Searches for all images paths and creates a job file

    Parameters
    ----------
    jobfile_path: Path to save job file to
    root_dir: the root project directory
    is_mutants: if True search the folder for individual line sub folders

    """
    # output_dir = root_dir / 'radiomics_output'
    # output_dir.mkdir(exist_ok=True)

    jobs_entries = []
    # get each file path
    for i, vol_path in enumerate(file_paths):

        rel_path_to_specimen_input = str(Path(vol_path).relative_to(jobs_file.parent / "rigids"))
        jobs_entries.append([rel_path_to_specimen_input,'to_run', '_', '_', '_'])

    jobs_df = pd.DataFrame.from_records(jobs_entries, columns=['job', 'status', 'host', 'start_time', 'end_time'])

    jobs_df.to_csv(jobs_file)
    return True


def denoise(images):
    ''''Lets just try out a denoiser'''

    #out_dir = _dir / "Patched_denoised"


    denoise = sitk.PatchBasedDenoisingImageFilter()

    for i, img in enumerate(images):


        img_arr = sitk.GetArrayFromImage(img)

        img_to_denoise = sitk.GetImageFromArray(img_arr)

        #cropped_arr = img_arr[100:300, 100:300, 100:300]



        #cropped_img.CopyInformation(img)
        logging.info("PIZZA TIME!")
        images[i] = denoise.Execute(img_to_denoise).CopyInformation(img)
    return images



def pyr_normaliser(_dir, _normaliser, scans_imgs, masks: list = None, fold: bool = False, ref_vol_path: Path = None, stage_for_ref: bool = False):
    # create a copy so orginal files aren't overwritten

    # Do the normalisation
    if isinstance(_normaliser, normalise.NonRegMaskNormalise):
        # checks if a ref mean has been calculated and then creates if missing
        if not _normaliser.reference_mean:
            #if you passed a non-normal label for reference
            if (ref_vol_path and stage_for_ref):
                ref_vol = common.LoadImage(ref_vol_path).img

                ref_mask_dir = ref_vol_path.parent.parent / "stage_labels"
                ref_mask_path = ref_mask_dir / os.path.basename(ref_vol_path)
                ref_mask = common.LoadImage(ref_mask_path).img


            elif (ref_vol_path):
                ref_vol = common.LoadImage(ref_vol_path).img
                ref_mask = _normaliser.gen_otsu_masks(ref_vol)

            else:
                # this one has multiple labels and volums but its singular to cut lines of code
                ref_vol, ref_mask = _normaliser.get_all_wt_vols_and_masks(_dir)

            _normaliser.add_reference(ref_vol, ref_mask)

        _normaliser.normalise(scans_imgs, masks, fold=fold, temp_dir=_dir)

    elif isinstance(_normaliser, normalise.IntensityHistogramMatch):
        if ref_vol_path:
            logging.info(f"Using {ref_vol_path} as the reference image")
            ref_vol = common.LoadImage(ref_vol_path).img
            _normaliser.normalise(scans_imgs, ref_vol)

        else:
            _normaliser.normalise(scans_imgs, scans_imgs[0])

    elif isinstance(_normaliser, normalise.IntensityN4Normalise):
        otsu_masks = _normaliser.gen_otsu_masks(scans_imgs)
        _normaliser.normalise(scans_imgs, otsu_masks)

    return scans_imgs


def pyr_calc_all_features(img, lab, name, labs_of_int, spherify=None):
    full_results = pd.Series([])
    #lab.CopyInformation(img)

    arr = sitk.GetArrayFromImage(lab)

    if spherify:  # can be used as a control - makes label a sphere:
        logging.info("Spherifying")
        s = ndimage.find_objects(arr)[-1]
        if spherify == 0:
            sphere_dir = Path(name).parent.parent / "colat_tumours"
            os.makedirs(sphere_dir, exist_ok=True)
            logging.info("Placing tumour as colateral control")
            lab = sitk.Flip(lab, [False, False, True])
            sphere_fname = sphere_dir / os.path.basename(name)
            sitk.WriteImage(lab, str(sphere_fname))

        elif spherify == 1:
            sphere_dir = Path(name).parent.parent / "spheres"
            os.makedirs(sphere_dir, exist_ok=True)
            logging.info("Spherifying in centre of tumour")
            midpoint = [np.round(np.mean([s[0].start, s[0].stop])) / 512,
                        np.round((np.mean([s[1].start, s[1].stop]))) / 512,
                        np.round((np.mean([s[2].start, s[2].stop]))) / 512]
            arr = rg.sphere(512, 10, midpoint, smoothing=True).astype(np.int_)
            mask = sitk.GetImageFromArray(arr)
            sphere_fname = sphere_dir / os.path.basename(name)
            sitk.WriteImage(mask, str(sphere_fname))

        else:
            sphere_dir = Path(name).parent.parent / "lateral_spheres"
            os.makedirs(sphere_dir, exist_ok=True)
            logging.info("Spherifying as colateral control")
            midpoint = [np.round(np.mean([s[0].start, s[0].stop])) / 512,
                        np.round((np.mean([s[1].start, s[1].stop]))) / 512,
                        np.round(482 - (np.mean([s[2].start, s[2].stop]))) / 512]
            arr = rg.sphere(512, 10, midpoint, smoothing=True).astype(np.int_)
            mask = sitk.GetImageFromArray(arr)
            sphere_fname = sphere_dir / os.path.basename(name)
            sitk.WriteImage(mask, str(sphere_fname))

    extractor = featureextractor.RadiomicsFeatureExtractor()
    extractor.enableAllImageTypes()
    extractor.enableAllFeatures()

    results_list =[]
    # TODO: reduce dimensionality?
    for i, org in enumerate(labs_of_int):
        # remove other labels

        arr_spec = np.where(arr == org, 1, 0)

        if np.count_nonzero(arr_spec) < 1000:
            print("null label")
            continue

        mask = sitk.GetImageFromArray(arr_spec)


        # make sure its in the same orientation as the image
        mask.CopyInformation(lab)


        result = extractor.execute(img, mask)

        features = pd.DataFrame.from_dict(result, orient='index',
                                          columns=[org]).transpose()

        #features = features.

        features = features.drop(columns=[col for col in features.columns if 'diagnostics' in col])
        #features = features.T.rename(columns={0: org})
        results_list.append(features)

    full_results = pd.concat(results_list, axis=0)

    return full_results


def run_radiomics(rad_dir, rigids, labels, name, labs_of_int,
                  norm_method, norm_label=None, spherify=None, ref_vol_path=None):
    """

    Parameters
    ----------
    norm_label : object
    """
    logging.info(common.git_log())
    signal.signal(signal.SIGINT, common.service_shutdown)
    mem_monitor = MonitorMemory(Path(rad_dir).absolute())

    features = pyr_calc_all_features(rigids, labels, name, labs_of_int, spherify=spherify)

    feature_dir = rad_dir / "features"

    os.makedirs(feature_dir, exist_ok=True)

    file_name = feature_dir / str(str(os.path.splitext(os.path.basename(name))[0]) + '.csv')

    features.to_csv(file_name)

    mem_monitor.stop()
    return True


def radiomics_job_runner(target_dir, labs_of_int=None,
                         norm_method=normalise.IntensityN4Normalise(),
                         norm_label=None, spherify=None,
                         ref_vol_path=None,
                         make_job_file: bool=False, fold: bool=False):
    '''
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
    # create files if they don't exist
    rad_dir = target_dir / 'radiomics_output'

    jobs_file_path = rad_dir / JOBFILE_NAME
    lock_file = jobs_file_path.with_suffix('.lock')
    lock = SoftFileLock(lock_file)

    if make_job_file:
        # extract the registrations if the job file doesn't exist and normalise
        if not os.path.exists(str(rad_dir)):
            os.makedirs(rad_dir, exist_ok=True)
            logging.info("Extracting Rigids")
            rigids = extract_registrations(target_dir)
            logging.info("Extracting Inverted Labels")
            labels = extract_registrations(target_dir, labs_of_int)

            if norm_label:
                logging.info("Extracting Stage labels")
                stage_labels = extract_registrations(target_dir, labs_of_int, norm_label=True)
            else:
                logging.info("Extracting Inverted Stats Masks")
                inv_stats_masks = extract_registrations(target_dir, labs_of_int, stats_mask=True)

        else: # good for debugging if normalisation stuffs up
            logging.info("loading rigids")
            rigids = [common.LoadImage(path).img for path in common.get_file_paths(str(rad_dir / "rigids"))]
            # labels = [common.LoadImage(path) for path in common.get_file_paths(str(rad_dir / "inverted_labels"))]
            logging.info("loading stats masks")
            inv_stats_masks = [common.LoadImage(path).img for path in common.get_file_paths(str(rad_dir / "stats_mask"))]
            stage_labels = [common.LoadImage(path).img for path in common.get_file_paths(str(rad_dir / "stage_labels"))]

        #logging.info("Denoising")

        #denoise(rigids)

        #logging.info("Writing Denoised Rigids")
        #rigid_paths = [common.LoadImage(path).img_path for path in common.get_file_paths(str(rad_dir / "rigids"))]
        # sort should be identical:
        # rigid_paths.sort(key=lambda x: os.path.basename(x))

        #for i, vol in enumerate(rigids):
        #    logging.info("Writing: {}".format(rigid_paths[i]))
        #    sitk.WriteImage(vol, rigid_paths[i], useCompression=True)



        # Normalisation should be here!!!!
        logging.info("Normalising Intensities")

        def prepare_norm(meth, rigids):
            if norm_label:
                rigids = pyr_normaliser(rad_dir, norm_method, scans_imgs=rigids, masks=stage_labels,
                                        ref_vol_path=ref_vol_path,
                                        stage_for_ref=True, fold=fold)
            else:
                if isinstance(meth, normalise.NonRegMaskNormalise):
                    logging.info("Normalising based on inverted stats masks")
                    rigids = pyr_normaliser(rad_dir, meth, scans_imgs=rigids, masks=inv_stats_masks)
                else:
                    rigids = pyr_normaliser(rad_dir, meth, scans_imgs=rigids, ref_vol_path=ref_vol_path)


        if isinstance(norm_method, list):
            for meth in norm_method:
                prepare_norm(meth, rigids)
        else:
            prepare_norm(norm_method, rigids)



        logging.info("Writing Normalised Rigids")
        rigid_paths = [common.LoadImage(path).img_path for path in common.get_file_paths(str(rad_dir / "rigids"))]
        # sort should be identical:
        #rigid_paths.sort(key=lambda x: os.path.basename(x))

        for i, vol in enumerate(rigids):
            logging.info("Writing: {}".format(rigid_paths[i]))
            sitk.WriteImage(vol, rigid_paths[i], useCompression=True)

        logging.info("Creating a job-file for radiomics")
        make_rad_jobs_file(jobs_file_path, rigid_paths)
        logging.info("Job_file_created")
        return True

    #df_jobs = pd.read_csv(jobs_file_path, index_col=0)

    # execute parallelisation:
    while True:
        try:
            with lock.acquire(timeout=60):

                df_jobs = pd.read_csv(jobs_file_path, index_col=0)
                # Get an unfinished job
                jobs_to_do = df_jobs[df_jobs['status'] == 'to_run']
                if len(jobs_to_do) < 1:
                    logging.info("No more jobs left on jobs list")

                    # error trap for processes that hung
                    logging.info("checking for hung jobs")
                    fin_jobs = df_jobs[df_jobs['status'] == 'complete']
                    running_jobs = df_jobs[df_jobs['status'] == 'running']
                    fin_indx = fin_jobs.index[-1]
                    fin_t = fin_jobs.at[fin_indx, 'start_time']
                    fin_time = datetime.strptime(fin_t, '%Y-%m-%d %H:%M:%S')
                    run_t = running_jobs['start_time']
                    run_times = [datetime.strptime(t, '%Y-%m-%d %H:%M:%S') < fin_time for t in run_t]
                    hung_jobs = running_jobs[run_times]

                    if len(hung_jobs) > 0:
                        logging.info("Hung jobs found - rerunning")
                        jobs_to_do = hung_jobs
                    else:
                        break
                indx = jobs_to_do.index[0]

                img_path = Path(rad_dir / 'rigids') / (jobs_to_do.at[indx, 'job'])
                lab_path = Path(rad_dir / 'inverted_labels') / (jobs_to_do.at[indx, 'job'])

                img = common.LoadImage(img_path)

                lab = common.LoadImage(lab_path)

                df_jobs.at[indx, 'status'] = 'running'
                df_jobs.at[indx, 'start_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                df_jobs.at[indx, 'host'] = socket.gethostname()
                df_jobs.to_csv(jobs_file_path)

        except Timeout:
            sys.exit('Timed out' + socket.gethostname())

        # try:
        try:
            logging.info(f'trying {img.img_path} and {lab_path}')
            run_radiomics(rad_dir, img.img, lab.img, img.img_path,
                          labs_of_int, norm_method, norm_label=norm_label, spherify=spherify)

        except Exception as e:
            if e.__class__.__name__ == 'KeyboardInterrupt':
                logging.info('terminating')
                sys.exit('Exiting')

            status = 'failed'
            print(e)
            logging.exception(e)


        else:
            status = 'complete'

        finally:
            with lock:
                df_jobs = pd.read_csv(jobs_file_path, index_col=0)
                df_jobs.at[indx, 'status'] = status
                df_jobs.at[indx, 'end_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                df_jobs.to_csv(jobs_file_path)

    logging.info('Exiting job_runner')
    return True
