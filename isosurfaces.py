__author__ = 'james'

import argparse
import SimpleITK as sitk
from os.path import join, dirname
import vtk


def generate_isosurfaces(input_nrrd, out_dir, organ_names):

    # Create smoothed labelmap
    im = sitk.ReadImage(input_nrrd)
    labelmap = join(dirname(input_nrrd), 'seg_smoothed.nrrd')
    smooth_and_fill_label_map(im, labelmap)

    # Create NRRD reader
    reader = vtk.vtkNrrdReader()
    reader.SetFileName(labelmap)

    for i, organ in enumerate(organ_names):

        file_name = '{}.vtk'.format(organ)
        vtk_out = join(out_dir, file_name)
        make_mesh(i+1, reader, vtk_out)


def smooth_and_fill_label_map(im, out_path):

    # Morphological opening
    opening = sitk.BinaryMorphologicalOpeningImageFilter()
    opening.SetKernelRadius([2, 2, 2])
    opening.SetKernelType(sitk.sitkBall)
    opened = opening.Execute(im > 0)

    # Fill holes
    filler = sitk.BinaryFillholeImageFilter()
    filler.FullyConnectedOn()
    filled = filler.Execute(opened)

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

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='vol', required=True, help="Input labelmap e.g. nrrd")
    parser.add_argument('-o', '--output_dir', dest='out_dir', required=True, help="Output directory")
    parser.add_argument('-n', '--names', dest='organ_names', required=True, help="CSV file of organ names")
    args = parser.parse_args()

    generate_isosurfaces(args.vol, args.out_dir, args.organ_names)
