#!/usr/bin/env python


import argparse
import SimpleITK as sitk
from os.path import join, dirname
import vtk
import csv

import common
from invert import BatchInvertMeshes


def generate_isosurfaces(inversion_config, input_nrrd, out_dir, organs_csv):
    """
    Cannot take compressed nrrds at the moment as VTK does not support them

    Parameters
    ----------
    inversion_config: str
        path to yaml file with directory names containing the elastix inversion parameter files
        eg:
    """

    # Create smoothed labelmap
    im = sitk.ReadImage(input_nrrd)
    labelmap = join(dirname(input_nrrd), 'seg_smoothed.tiff')
    smooth_and_fill_label_map(im, labelmap)

    organ_names = []
    with open(organs_csv, 'r') as on:
            reader = csv.reader(on)
            for row in reader:
                try:
                    organ_names.append(row[0])
                except IndexError:
                    pass

    reader = vtk.vtkTIFFReader()
    reader.SetFileName(labelmap)

    mesh_dir = join(out_dir, 'meshes')
    common.mkdir_if_not_exists(mesh_dir)

    elastix_points_dir = join(out_dir, 'inverted_meshes_elx_points')
    common.mkdir_if_not_exists(elastix_points_dir)

    mesh_filenames = []
    for i, organ in enumerate(organ_names):

        file_name = '{}.vtk'.format(organ)
        mesh_filenames.append(file_name)
        vtk_out = join(mesh_dir, file_name)
        make_mesh(i+1, reader, vtk_out)
        BatchInvertMeshes(inversion_config, vtk_out, elastix_points_dir)


# def vtk_to_elastix_points(vtk_file, outfile):
#     """
#     Converts vtk points from a polygon file to a format that easltix can use
#     """
#
#     def points(l, n):
#         """
#         Yield successive n-sized chunks from l.
#         """
#         for i in xrange(0, len(l), n):
#             yield l[i:i+n]
#
#     point_lines = []
#     numpoints = 0
#
#     with open(vtk_file, 'r') as fin, open(outfile, 'w') as fout:
#         start_reading = False
#         for line in fin.readlines():
#             if start_reading:
#                 l = line.split()
#                 if line.startswith('POLYGONS'):
#                     break
#                 for p in points(l, 3):
#                     numpoints += 1
#                     point_lines.append(" ".join(p) + "\n")
#             elif line.startswith('POINTS'):
#                 start_reading = True
#
#         fout.write("{}\n".format(numpoints))
#         fout.writelines(point_lines)


def smooth_and_fill_label_map(im, out_path):

    # Morphological opening
    # opening = sitk.BinaryMorphologicalOpeningImageFilter()
    # opening.SetKernelRadius([1, 1, 1])
    # opening.SetKernelType(sitk.sitkBall)
    # opened = opening.Execute(im > 0)

    # Fill holes
    filler = sitk.BinaryFillholeImageFilter()
    filler.FullyConnectedOn()
    filled = filler.Execute(im)

    # Merge the mask with the original image
    masking = sitk.MaskImageFilter()
    im = masking.Execute(im, filled)

    # Save, uncompressed (!)
    sitk.WriteImage(im, out_path, False)


def make_mesh(label_num, reader, vtk_out):

    # Marching cubes
    mcsurf = vtk.vtkDiscreteMarchingCubes()
    mcsurf.SetInputConnection(reader.GetOutputPort())
    mcsurf.SetValue(0, label_num)
    mcsurf.ComputeScalarsOff()
    mcsurf.ComputeGradientsOff()

    # Smoothing
    smoothing = vtk.vtkSmoothPolyDataFilter()
    smoothing.SetInputConnection(mcsurf.GetOutputPort())
    smoothing.SetRelaxationFactor(0.1)
    smoothing.SetNumberOfIterations(100)

    # Decimate
    decimate = vtk.vtkDecimatePro()
    decimate.SetInputConnection(smoothing.GetOutputPort())
    decimate.SetTargetReduction(.9)  # 90% reduction

    # Write VTK
    vtk_writer = vtk.vtkPolyDataWriter()
    vtk_writer.SetInputConnection(decimate.GetOutputPort())
    vtk_writer.SetFileName(vtk_out)
    vtk_writer.Write()


def get_label_range(im):

    arr = sitk.GetArrayFromImage(im)
    min_ = arr.min()
    max_ = arr.max()
    return range(min_ + 1, max_)


def save_vtk(vtk_data, output_path):

    # Write shape as VTK
    vtk_writer = vtk.vtkPolyDataWriter()
    vtk_writer.SetFileName(output_path)
    vtk_writer.SetInputData(vtk_data)
    vtk_writer.Write()

if __name__ == '__main__':
    import sys
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', dest='config', required=True, help="Input labelmap e.g. nrrd")
    parser.add_argument('-i', '--invertable', dest='vol', required=True, help="Output directory")
    parser.add_argument('-o', '--out', dest='out_dir', required=True, help="CSV file of organ names")
    parser.add_argument('-n', '--names', dest='organ_names', required=True, help="yaml cofig with inversion order")

    args, _ = parser.parse_known_args()



    generate_isosurfaces(args.config, args.vol, args.out_dir, args.organ_names)
