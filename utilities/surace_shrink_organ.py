#!/usr/bin/env python

import vtk
import subprocess

import SimpleITK as sitk

NORMAL_ORGAN_MESH_PATH = 'organ.vtk'
SHRUNK_ORGAN_MESH_PATH = 'shrunk_organ.vtk'

ELASTIX_PATH = 'elxParamTemp.txt'

PARAM = """
(AutomaticScalesEstimation  "true")
(NumberOfHistogramBins  32)
(MovingImageDimension  3)
(DefaultPixelValue  0)
(UseDirectionCosines  "true")
(FixedImageDimension  3)
(NewSamplesEveryIteration  "true")
(FixedImagePyramid  "FixedSmoothingImagePyramid")
(FinalBSplineInterpolationOrder  3)
(Resampler  "DefaultResampler")
(WriteResultImage  "true")
(CompressResultImage  "true")
(ResultImageFormat  "nrrd")
(ImageSampler  "Random")
(WriteTransformParametersEachIteration  "false")
(FixedInternalImagePixelType  "short")
(MovingImagePyramid  "MovingSmoothingImagePyramid")
(Interpolator  "BSplineInterpolator")
(ResultImagePixelType  "short")
(AutomaticTransformInitialization  "true")
(MovingInternalImagePixelType  "short")
(HowToCombineTransforms  "Compose")
(ResampleInterpolator  "FinalBSplineInterpolator")
(AutomaticParameterEstimation  "true")
(Optimizer  "AdaptiveStochasticGradientDescent")
(FinalGridSpacingInVoxels  16)
(MaximumNumberOfIterations  150)
(ASGDParameterEstimationMethod  "DisplacementDistribution")
(Metric "AdvancedMattesMutualInformation" "CorrespondingPointsEuclideanDistanceMetric")
(Metric0Weight 0.0)
(Metric1Weight 1.0)
(Transform  "BSplineTransform")
(MaximumStepLength  2.0)
(NumberOfSpatialSamples  3000)
(Registration  "MultiMetricMultiResolutionRegistration")
(UseAdaptiveStepSizes  "true")
(NumberOfResolutions  1)
"""

def labelmap_to_isosurfaces(vol, img, reg_out_dir):
    """

    :param vol:
    :param outdir:
    :param decimation:
    :param smoothiong:
    :return:
    """
    reader = vtk.vtkTIFFReader()
    reader.SetFileName(vol)

    mesh = make_mesh(1, reader)
    shrink_mesh(mesh)

    fixed_points = 'shrink_organ_points.elx'
    moving_points = 'normal_organ_points.elx'

    vtk_to_elastix_pounts(NORMAL_ORGAN_MESH_PATH, moving_points)
    vtk_to_elastix_pounts(SHRUNK_ORGAN_MESH_PATH, fixed_points)

    warp_using_points(img, fixed_points, moving_points, reg_out_dir)

def warp_using_points(img, fixed_points, moving_points, out_dir):


        with open(ELASTIX_PATH, 'w') as fh:
            fh.write(PARAM)

        cmd = ['elastix',
               '-f', img,
               '-fp', fixed_points,
               '-m', img,
               '-mp', moving_points,
               '-out', out_dir,
               '-p', ELASTIX_PATH
              ]

        subprocess.check_output(cmd)


def make_mesh(label_num, reader):

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

    # connectivity filter
    cfsurf = vtk.vtkPolyDataConnectivityFilter()
    cfsurf.SetExtractionModeToLargestRegion()
    cfsurf.SetInputConnection(mcsurf.GetOutputPort())
    cfsurf.InitializeSpecifiedRegionList()

    cfsurf.Update()
    cleanfilter = vtk.vtkCleanPolyData()
    cleanfilter.SetInputConnection(cfsurf.GetOutputPort())
    cleanfilter.Update()

    mcsurf.SetValue(0,100)

    writer = vtk.vtkPolyDataWriter()
    writer.SetInputConnection(cfsurf.GetOutputPort())
    writer.SetFileName(NORMAL_ORGAN_MESH_PATH)
    writer.Write()

    # At the moment I'm writing the output from vtkPolyDataConnectivityFilter as I don't know how to get the points



    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(NORMAL_ORGAN_MESH_PATH)
    reader.Update()

    poly = reader.GetOutput()

    return poly

def shrink_mesh(polydata):
    SCALE = 0.9

    com_filter = vtk.vtkCenterOfMass()
    com_filter.SetInputData(polydata)
    com_filter.SetUseScalarsAsWeights(False)
    com_filter.Update()

    cx, cy, cz = com_filter.GetCenter()

    #Do transform

    tform = vtk.vtkTransform()
    tform.Translate(cx, cy, cz)
    tform.Scale(SCALE, SCALE, SCALE)
    tform.Translate(-cx, -cy, -cz)

    tfilter = vtk.vtkTransformPolyDataFilter()
    tfilter.SetInputData(polydata)
    tfilter.SetTransform(tform)
    tfilter.Update()

    writer = vtk.vtkPolyDataWriter()
    writer.SetInputConnection(tfilter.GetOutputPort())
    writer.SetFileName(SHRUNK_ORGAN_MESH_PATH)
    writer.Write()


def vtk_to_elastix_pounts(in_, out):

    reading = False
    xyz_points = []
    with open(in_, 'r') as fin, open(out, 'w') as fout:
        for line in fin.readlines():
            l = line.split()
            if not reading:
                try:
                    float(l[0])
                except ValueError:
                    continue
                else:
                    reading = True

            if reading:
                try:
                    if l[0].lower().startswith('polygons'):
                         break
                except IndexError:
                    print 'something'

                for p in points(l, 3):
                    xyz_points.append(p)
                    #fout.write(" ".join(p))
                    #fout.write('\n')
        fout.write('index\n')
        fout.write(str(len(xyz_points)) + '\n')
        for xyz in xyz_points:
            fout.write(" ".join(xyz) + '\n')


def points(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


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
    img = sys.argv[2]
    reg_out_dir = sys.argv[3]
    labelmap_to_isosurfaces(labelMap, img, reg_out_dir)