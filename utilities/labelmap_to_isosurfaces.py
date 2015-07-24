from scipy.weave.blitz_spec import array_converter

__author__ = 'james'

import vtk
import argparse
import SimpleITK as sitk
import os
import numpy as np

def labelmap_to_isosurfaces(vol, outdir, decimation, smoothing):
    """

    :param vol:
    :param outdir:
    :param decimation:
    :param smoothiong:
    :return:
    """
    print vol
    reader = vtk.vtkTIFFReader()
    reader.SetFileName(vol)



    mesh = make_mesh(1, reader, outdir)
    return mesh


def make_mesh(label_num, reader, outdir):

    imgthresh = vtk.vtkImageThreshold()
    imgthresh.SetInputConnection(reader.GetOutputPort())
    imgthresh.ThresholdBetween(label_num, label_num)
    imgthresh.ReplaceInOn()
    imgthresh.SetInValue(100)
    imgthresh.SetOutValue(0)
    imgthresh.Update()

    mcsurf = vtk.vtkMarchingCubes()
    mcsurf.SetInputConnection(imgthresh.GetOutputPort())
    mcsurf.ComputeScalarsOff()
    mcsurf.ComputeGradientsOff()

    # desurf = vtk.vtkDecimatePro()
    # desurf.SetInputConnection(mcsurf.GetOutputPort())
    # #desurf.SetTargetReduction(0.1)

    # remember don't use the stripper

    # connectivity filter
    cfsurf = vtk.vtkPolyDataConnectivityFilter()
    cfsurf.SetExtractionModeToLargestRegion()
    cfsurf.SetInputConnection(mcsurf.GetOutputPort())
    cfsurf.InitializeSpecifiedRegionList()



    cfsurf.Update()
    cleanfilter = vtk.vtkCleanPolyData()
    cleanfilter.SetInputConnection(cfsurf.GetOutputPort())
    cleanfilter.Update()
    p = cleanfilter.GetOutput().GetNumberOfPoints()
    print 'olkdjkcflojkds', p

    com_filter = vtk.vtkCenterOfMass()
    com_filter.SetInputData(cfsurf.GetOutput())
    com_filter.SetUseScalarsAsWeights(False)
    com_filter.Update()


    com = com_filter.GetCenter()

    # print centroid

    #Do the shrinktf2 = vtk.vtkTransformPolyDataFilter()
    shrinkFilter = vtk.vtkShrinkPolyData()
    shrinkFilter.SetShrinkFactor(0.5)
    shrinkFilter.SetInputConnection(cfsurf.GetOutputPort())
    shrinkFilter.Update()


    mcsurf.SetValue(0,100)

    writer = vtk.vtkPolyDataWriter()
    writer.SetInputConnection(cfsurf.GetOutputPort())
    writer.SetFileName(os.path.join(outdir, str(label_num) + ".vtk"))
    writer.Write()

    writer.SetInputConnection(shrinkFilter.GetOutputPort())
    writer.SetFileName(os.path.join(outdir, str(label_num) + "shrunk.vtk"))
    writer.Write()


def array_to_vtkmatrix4x4(scipy_array):
            mat = vtk.vtkMatrix4x4()
            for i in range(0,4):
                for j in range(0,4):
                    mat.SetElement(i,j, scipy_array[i][j])
            return mat

def getlabelrange(vol):
    """
    Get the min and max label value from a label map and return the sequence from min to max
    :return:
    """
    img = sitk.ReadImage(vol)
    arr = sitk.GetArrayFromImage(img)
    min_ = arr.min()
    max_ = arr.max()
    return range(min_, max_ + 1)


if __name__ == '__main__':
    import sys
    labelMap = sys.argv[1]
    outfir = sys.argv[2]
    labelmap_to_isosurfaces(labelMap, outfir, 2, 2)