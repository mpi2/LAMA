import SimpleITK as sitk
import os
import numpy as np


def GetFilePaths(folder, extension_tuple=('.nrrd', '.tiff', '.tif', '.nii', '.bmp', 'jpg '), pattern=None):
    """
    Test whether input is a folder or a file. If a file or list, return it.
    If a dir, return all images within that directory.
    Optionally test for a pattern to sarch for in the filenames
    """
    if not os.path.isdir(folder):
        if isinstance(folder, basestring):
            return [folder]
        else:
            return folder
    else:
        paths = []
        for root, _, files in os.walk(folder):
            for filename in files:
                if filename.endswith(extension_tuple):
                    if pattern:
                        if pattern and pattern not in filename:
                            continue
                    #paths.append(os.path.abspath(os.path.join(root, filename))) #Broken on shared drive
                    paths.append(os.path.join(root, filename))
        return paths


def ItkImageGenerator(folder, extension_tuple=('.nrrd', '.tiff', '.nii', '.bmp')):
    '''
    Given a foler path and a list of possible file extensions
    Return all paths to found files including in subdirectoris
    '''
    for im_path in GetFilePaths(folder):
        itk_img = sitk.ReadImage(im_path)
        yield itk_img


def ImagesToNumpyArrays(folder, extension_tuple=('.nrrd', '.tiff', '.nii', '.bmp')):
    '''
    Given a foler path and a list of possible file extensions
    Return all paths to found files including in subdirectoris
    '''
    for im_path in GetFilePaths(folder):
        itk_img = sitk.ReadImage(im_path)
        yield sitk.GetArrayFromImage(itk_img)


def Average(img_dirOrList, search_subdirs=True):
    '''
    Create an average volume from multiple volumes
    @return: sitk Image
    '''
    images = []
    if isinstance(img_dirOrList, basestring):
        images = GetFilePaths(img_dirOrList)
    else:
        images = img_dirOrList

    #sum all images together
    first_img = sitk.ReadImage(images[0])
    summed = sitk.GetArrayFromImage(first_img)
    for image in images[1:]:#Ommit the first as we have that already
        itk_img = sitk.ReadImage(image)
        np_array = sitk.GetArrayFromImage(itk_img)
        try:
            summed = summed + np_array
        except ValueError as e:
            print "Numpy can't average this volume {0}".format(image)

    #Now make average. Do it in numpy as I know how
    avg_array = summed/len(images)
    avg_img = sitk.GetImageFromArray(avg_array)
    return avg_img


# def CastTo8BitScale255(file_or_dir, outdir):
#     if os.path.isdir(file_or_dir):
#         filelist = GetFilePaths(file_or_dir)
#     else:
#         filelist = [file_or_dir]
#     for im_path in filelist:
#         img = sitk.ReadImage(im_path)
#         scaled_img = sitk.RescaleIntensity(img, 0.0, 255.00)
#         bit8_img = sitk.Cast(scaled_img, sitk.sitkUInt8)
#         AppendToFileBaseName(im_path, "_scaled" )
#         sitk.WriteImage(bit8_img, AppendToFileBaseName(im_path, "_scaled" ) )#Change this to yield images

def rescale_cast8bit(volumes, outpath):
    if isinstance(volumes, basestring):
        if os.path.isdir(volumes):
            filelist = GetFilePaths(volumes)
        else:
            filelist = [volumes]
    else:
        filelist = volumes # already a list
    for fn in filelist:
        img = sitk.ReadImage(fn)
        arr = sitk.GetArrayFromImage(img)
        min = arr.min()
        max = arr.max()
        if min > 0:
            print "min of {fn} not less than 0, maybe doesn't need scaling".format(**locals())
        scaled = sitk.RescaleIntensity(img, 0.0, float(max+abs(min)))
        scaled_cast = sitk.Cast(scaled, sitk.sitkUInt8)
        basename = os.path.splitext(os.path.basename(fn))[0]
        outfn = os.path.join(outpath, basename) + 'scaled_8bit.nrrd'
        sitk.WriteImage(scaled_cast, outfn )


def AppendToFileBaseName(path, addition):
    base, ext = os.path.splitext(os.path.basename(path))
    path = os.path.split(path)[0]
    new_base = base + addition
    return os.path.join(path, new_base + ext)


def PrependToFileBasename(path):
    pass



