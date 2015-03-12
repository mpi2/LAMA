#!/usr/bin/env python

"""
Generate a statistics data file that can be used as input for vpv.py
"""
from simplejson import ordered_dict

import numpy as np
import SimpleITK as sitk
import argparse
import cPickle as pickle
from reg_pipeline.process import harwellimglib as hil
import sys
import os
from multiprocessing import Pool
import pprint
from collections import OrderedDict


def make_cubes(process_args):
    """
    chop a volume/dewformation file into x3 sized chunks, return a list with top_left pos, mean magnitude, and filename
    :return:mean_cubes: list, [[(z,y,x), mean_vect_magnitude(float),filename],.....]
    """
    global chunksize
    data_type, rawdata_file = process_args

    filename = os.path.basename(rawdata_file)
    img = sitk.ReadImage(rawdata_file)
    if data_type == 'intensities':
        print('normalizing')
        img = sitk.Normalize(img)


    array = sitk.GetArrayFromImage(img)
    zdim, ydim, xdim = array.shape[0:3]  # For vectors there's and extra dimension so can't jsut unpack

    mean_cubes = []

    for z in range(0, zdim - chunksize, chunksize):
        for y in range(0, ydim - chunksize, chunksize):
            for x in range(0, xdim - chunksize, chunksize):
                cubedata = [(z, y, x)]  # start coordinates of cube
                if data_type == 'deformation_vectors':
                    mean_magnitude = cube_vect_magnitude_mean(array[z:z+chunksize, y:y+chunksize, x:x+chunksize])
                else:
                    mean_magnitude = cube_jacobian_mean(array[z:z+chunksize, y:y+chunksize, x:x+chunksize])
                cubedata.append(mean_magnitude)
                cubedata.append(filename)
                mean_cubes.append(cubedata)

    print(filename, 'done')
    return mean_cubes  #[(x,y,z), [data,..], [x,y..]


def cube_vect_magnitude_mean(cube_of_vectors):
    """
    For a cube of deformation vectors, get the mean magnitude
    :param cube_of_vectors:
    :return: mean magnitude
    """
    vectors = []
    #Append each vector from the cube to a list
    for z in cube_of_vectors:
        for y in z:
            for vec in y:
                vectors.append(vec)
    #Get the mean vector. Then get the magnitude of it using np.linalg.norm
    return np.linalg.norm(np.mean(vectors, axis=0))


def cube_jacobian_mean(cube_of_jac_scalars):
    """
    For a cube of jacobian scalars, get the mean value
    :param cube_of_jac_scalars:
    :return:mean jacobian
    """
    return np.mean(cube_of_jac_scalars)



def group_position_means(all_files_data):
    """
    Ggroup together mean magnitudes from each file by corresponding position.
    :param all_files_data: Data for all files at one cube position
    :return:pos_means, list
    """
    pos = [all_files_data[0][0]]  #=(x,y,z)
    all_means = []

    #add the each files cube average
    for deform_file_data in all_files_data:
        all_means.append(deform_file_data[1])

    pos.append(all_means)
    return pos


def run(args):
    """
    :param jacobians:List of jacobian files
    :param deforms: List of defomation files
    :param outfile: Where to put the output
    :param threads: Max num of threads to use
    :param deform_file_pattern: defomation file regex finder in case these files are mixed in with others
    """

    if os.path.isdir(args.outfile):
        sys.exit("Supply an output file path not a directory")
    global chunksize

    if args.deforms:  # Get deformation fields
        print('processing deformation vectors')
        def_files = hil.GetFilePaths(os.path.abspath(args.deforms), pattern=args.deform_file_pattern)
        if len(def_files) < 1:
            sys.exit("No files found in the directory: {}".format(args.deforms))
        average_deformation(def_files, args.threads, 'deformation_vectors', args.outfile + '_def_vec')

    if args.jacobians:  # Get jacobians
        print('processing jacobians')
        jac_files = hil.GetFilePaths(os.path.abspath(args.jacobians), pattern=args.jac_file_pattern)
        if len(jac_files) < 1:
            sys.exit("No files found in the directory: {}".format(args.jacobians))
        average_deformation(jac_files, args.threads, 'jacobians', args.outfile + '_jac')

    if args.registered_vols:  # Get intensities
        print('processing intesity differences')
        jac_files = hil.GetFilePaths(os.path.abspath(args.registered_vols))
        if len(jac_files) < 1:
            sys.exit("No files found in the directory: {}".format(args.registered_vols))
        average_deformation(jac_files, args.threads, 'intensities', args.outfile + '_int_diff')


def average_deformation(files, threads, datatype, outfile):
    """.. function:: format_exception(etype, value, tb[, limit=None])

   Format the exception with a traceback.

   :param etype: exception type
   :param value: exception value
   :param tb: traceback object
   :param limit: maximum number of stack frames to show
   :type limit: integer or None
   :rtype: list of strings
    """
    #get the dimemsions of the deform files
    im = sitk.ReadImage(files[0])
    imarray = sitk.GetArrayFromImage(im)
    deform_dimensions = imarray.shape
    pool = Pool(processes=threads)

    # Multiprocess.map only accepts one iterable. so bundle up the filenames with the datatype and intensity ref'
    proc_args = [(datatype, volpath) for volpath in files]

    file_means = pool.map(make_cubes, proc_args)  # get list of lists containing each file's mean cube data.
    pos_means = map(group_position_means, zip(*file_means))   # Now group together the means from each position
    #Now add to a dict
    final_average_data = OrderedDict()
    for pm in pos_means:
        final_average_data[pm[0]] = pm[1]

    # fh = open(outfile + ".readable", 'wb')
    # #fh.write("location, cube means,\n")
    # pprint.pprint(cube_means, fh)
    # fh.close()
    # pickle.dump(cube_means, open(outfile, 'wb'), -1)


    #final_average_data = []
    # for cube in cube_means:
    #     final_average_data.append(cube)  # This seems useless. Just assigning one list to another ?

    deformation_output = OrderedDict({
        'data_type': datatype,
        'chunksize':  chunksize,
        'deform_dimensions': deform_dimensions,
        'vols': files,
        'data': final_average_data,
    })

    fh = open(outfile + "_readable", 'wb')
    #fh.write("location, cube means,\n")
    pprint.pprint(deformation_output, fh)
    fh.close()
    pickle.dump(deformation_output, open(outfile, 'wb'))


if __name__ == '__main__':

    # mpl = multiprocessing.log_to_stderr()
    # mpl.setLevel(multiprocessing.SUBDEBUG)

    parser = argparse.ArgumentParser("messageS")
    parser.add_argument('-c', dest='cubesize', type=int, help='Size of the sub array', required=True)
    parser.add_argument('-d', dest='deforms', help='Folder containing deformation field files ', default=None)
    parser.add_argument('-j', dest='jacobians', help='Folder containing jacobian field files ', default=None)
    parser.add_argument('-r', dest='registered_vols', help='Folder containing registered vols, for intensity difference'
                                                           ' analysis', default=None)
    parser.add_argument('-o', dest='outfile', help='File to store pickle file of means/stdvs of vectors,intensities etc', required=True)
    parser.add_argument('-t', dest='threads', type=int, help='How many threads to use', default=4)
    parser.add_argument('-dp', dest='deform_file_pattern', help='String that is contained in the deform file names',
                        default='deformationField')
    parser.add_argument('-jp', dest='jac_file_pattern', help='String that is contained in the jacobian file names',
                        default='spatialJacobian')

    args = parser.parse_args()
    if not args.deforms and not args.jacobians and not args.registered_vols:
        sys.exit('You need to supply deformations, jacobians, or registered volumes (-d, -j, -r)')

    #cProfile.run(run(args.deforms, args.cubesize, args.outfile))
    global chunksize
    chunksize = args.cubesize
    run(args)
