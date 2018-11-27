import numpy as np
import re
import math
import os
import sys
import subprocess
import shutil
import fileinput

try :
  import SimpleITK as sitk
  sitkLoaded = True;
except ImportError :
  sitkLoaded = False
  pass

#=========================================================================
#
#  Copyright Roy van Pelt - Eindhoven University of Technology (2013)
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#
#  Population-based registration of spheres example
#
#  Based on:
#  W. Van Hecke, J. Sijbers, E. D'Agostino, F. Maes, S. De Backer,
#  E. Vandervliet, P. M. Parizel, and A. Leemans
#  On the construction of an inter-subject diffusion tensor
#  magnetic resonance atlas of the healthy human brain.
#  NeuroImage, vol. 43, no. 1, pp. 69-80, Oct. 2008.
#
#  Amongst others also found in:
#  D. Seghers, E. D'Agostino, F. Maes, D. Vandermeulen, and P. Suetens.
#  Construction of a brain template from MR images using
#  state-of-the-art registration and segmentation techniques.
#  In Proc Intl Conf Med Image Comput Comp Ass Interv, pp 696-703. 2004.
#
#  Note:
#  To compute the transformation from the selected subject to the
#  atlas space, we take the inverse of the mean transformation.
#  Van Hecke et al. use the mean of all inverse transformations,
#  which likely leads to erroneous results.
#
#  Acknowledgements:
#  I'd like to thanks Marius Staring for the support on Elastix,
#  and atlas registration in general.
#
#  Using:
#  Python    v2.7.4
#  Elastix   v4.600 (http://elastix.isi.uu.nl/)
#  SimpleITK v0.6.1 (http://www.simpleitk.org/)
#
#=========================================================================

print ("==============================================")
print ("= Spheres - Population-based registration    =")
print ("= Roy van Pelt (2013)                        =")
print ("= Eindhoven University of Technology         =")
print ("==============================================\n")

# ------------------------------------------------------
# Current working directory
# ------------------------------------------------------
currentDir        = os.getcwd()

# For cygwin:
m = re.match( "/cygdrive/[a-z]/", currentDir )
if m :
  find =    m.group(0);
  replace = m.group(0).replace( "/cygdrive/", "" ).replace( "/", ":/" ).capitalize();
  currentDir = currentDir.replace( find, replace )

# ------------------------------------------------------
# User variables
# ------------------------------------------------------
outputDir         = os.path.join(currentDir, "output/")
parameterDir      = os.path.join(currentDir, "parameters/")
dataDir           = "spheres"

parameterFile0    = os.path.join(parameterDir, "parameters_sphere_affine.txt")
parameterFile1    = os.path.join(parameterDir, "parameters_sphere_nonrigid.txt")
parameterFileInv  = os.path.join(parameterDir, "parameters_sphere_inverse.txt")

scalarFile        = "sphere"
centersFile       = "sphere_centers"

dataSetCount      = 6

# ------------------------------------------------------
# Initialize variables
# ------------------------------------------------------
transformLabel0   = "TransformParameters.0.txt"
transformLabel1   = "TransformParameters.1.txt"
transformLabelIdt = "TransformParameters.idt.txt"

resultLabel0      = "result.0.mhd"
resultLabel1      = "result.1.mhd"
resultLabel       = "result.mhd"

scalarVolumes = []

for i in range(1,dataSetCount+1):
  scalarIndex = str(i).rjust(3, '0')
  scalarPath  = os.path.join(currentDir, "%s/%s_%s.mhd" % (dataDir, scalarFile, scalarIndex))
  scalarVolumes.append(scalarPath)

# ------------------------------------------------------
# Check for input data
# ------------------------------------------------------
for i in range(0,dataSetCount):
  if not os.path.exists( scalarVolumes[i] ) :
    print((">> Error: input volume " + scalarVolumes[i] + " cannot be found."))
    raise SystemExit

# ------------------------------------------------------
# Function definitions
# ------------------------------------------------------
def createTransformIdentity(fileName, dimensions = [], spacing = [], origin = [], transform = []):

  fileIdt = open(transformFileIdt, 'wb')

  fileIdt.write('(Transform \"TranslationTransform\")\n')
  fileIdt.write('(NumberOfParameters 3)\n')
  fileIdt.write('(TransformParameters 0.0 0.0 0.0)\n')
  fileIdt.write('(InitialTransformParametersFileName \"NoInitialTransform\")\n')
  fileIdt.write('(HowToCombineTransforms \"Compose\")\n\n')
  fileIdt.write('// Image Specific\n')
  fileIdt.write('(FixedImageDimension 3)\n')
  fileIdt.write('(MovingImageDimension 3)\n')
  fileIdt.write('(FixedInternalImagePixelType \"float\")\n')
  fileIdt.write('(MovingInternalImagePixelType \"float\")\n')
  fileIdt.write('(Size '      + str(dimensions[0])   + ' ' + str(dimensions[1])  + ' ' + str(dimensions[2]) + ')\n')
  fileIdt.write('(Index 0 0 0)\n')
  fileIdt.write('(Spacing '   + str(spacing[0])      + ' ' + str(spacing[1])      + ' ' + str(spacing[2]) + ')\n')
  fileIdt.write('(Origin '    + str(origin[0])       + ' ' + str(origin[1])       + ' ' + str(origin[2])  + ')\n')
  fileIdt.write('(Direction ' + str(transform[0][0]) + ' ' + str(transform[0][1]) + ' ' + str(transform[0][2]) + ' ' \
                              + str(transform[1][0]) + ' ' + str(transform[1][1]) + ' ' + str(transform[1][2]) + ' ' \
                              + str(transform[2][0]) + ' ' + str(transform[2][1]) + ' ' + str(transform[2][2]) + ')\n')
  fileIdt.write('(UseDirectionCosines \"true\")\n\n')
  fileIdt.write('// ResampleInterpolator specific\n')
  fileIdt.write('(ResampleInterpolator \"FinalBSplineInterpolator\")\n')
  fileIdt.write('(FinalBSplineInterpolationOrder 3)\n\n')
  fileIdt.write('// Resampler specific\n')
  fileIdt.write('(Resampler \"DefaultResampler\")\n')
  fileIdt.write('(DefaultPixelValue 0.000000)\n')
  fileIdt.write('(ResultImageFormat \"mhd\")\n')
  fileIdt.write('(ResultImagePixelType \"float\")\n')
  fileIdt.write('(CompressResultImage \"false\")\n')

  fileIdt.close()

def createTransformAverage(fileName, transforms, transformCount, transformWeights = [], dimensions = [], spacing = [], origin = [], transform = [] ):
  fileAvg = open(transformFileAvg, 'wb')

  fileAvg.write('(Transform \"WeightedCombinationTransform\")\n')
  fileAvg.write('(NumberOfParameters ' +str(transformCount)+')\n')
  fileAvg.write('(TransformParameters '+transformWeights+')\n')
  fileAvg.write('(InitialTransformParametersFileName \"NoInitialTransform\")\n')
  fileAvg.write('(NormalizeCombinationWeights \"true\")\n')
  fileAvg.write('(SubTransforms '      +transforms+')\n')
  fileAvg.write('(HowToCombineTransforms \"Compose\")\n\n')
  fileAvg.write('// Image Specific\n')
  fileAvg.write('(FixedImageDimension 3)\n')
  fileAvg.write('(MovingImageDimension 3)\n')
  fileAvg.write('(FixedInternalImagePixelType \"float\")\n')
  fileAvg.write('(MovingInternalImagePixelType \"float\")\n')
  fileAvg.write('(Size '      + str(dimensions[0])   + ' ' + str(dimensions[1])  + ' ' + str(dimensions[2]) + ')\n')
  fileAvg.write('(Index 0 0 0)\n')
  fileAvg.write('(Spacing '   + str(spacing[0])      + ' ' + str(spacing[1])      + ' ' + str(spacing[2]) + ')\n')
  fileAvg.write('(Origin '    + str(origin[0])       + ' ' + str(origin[1])       + ' ' + str(origin[2])  + ')\n')
  fileAvg.write('(Direction ' + str(transform[0][0]) + ' ' + str(transform[0][1]) + ' ' + str(transform[0][2]) + ' ' \
                              + str(transform[1][0]) + ' ' + str(transform[1][1]) + ' ' + str(transform[1][2]) + ' ' \
                              + str(transform[2][0]) + ' ' + str(transform[2][1]) + ' ' + str(transform[2][2]) + ')\n')
  fileAvg.write('(UseDirectionCosines \"true\")\n\n')
  fileAvg.write('// ResampleInterpolator specific\n')
  fileAvg.write('(ResampleInterpolator \"FinalBSplineInterpolator\")\n')
  fileAvg.write('(FinalBSplineInterpolationOrder 3)\n\n')
  fileAvg.write('// Resampler specific\n')
  fileAvg.write('(Resampler \"DefaultResampler\")\n')
  fileAvg.write('(DefaultPixelValue 0.000000)\n')
  fileAvg.write('(ResultImageFormat \"mhd\")\n')
  fileAvg.write('(ResultImagePixelType \"float\")\n')
  fileAvg.write('(CompressResultImage \"false\")\n')

  fileAvg.close()

# ------------------------------------------------------
# Compute suitable image characteristics for
# the output volume in the average domain
# ------------------------------------------------------
outputDimensions = [1,1,1]
outputSpacing    = [1,1,1]
outputOffset     = [0.0,0.0,0.0]
outputTransform  = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]

dims    = [[]*3 for i in range(3)]
spacing = [[]*3 for i in range(3)]
offset  = [[]*3 for i in range(3)]
axisX   = [[]*3 for i in range(3)]
axisY   = [[]*3 for i in range(3)]

for subject in range(1,dataSetCount+1):
  index    = str(subject).rjust(3, '0')
  volume   = os.path.join(currentDir, "%s/%s_%s.mhd" % (dataDir,scalarFile,index))

  for line in fileinput.input(volume):
    if (line.count('DimSize')>0) :
      dims[0].append(int(line.split()[2]))
      dims[1].append(int(line.split()[3]))
      dims[2].append(int(line.split()[4]))
    if (line.count('ElementSpacing')>0) :
      spacing[0].append(float(line.split()[2]))
      spacing[1].append(float(line.split()[3]))
      spacing[2].append(float(line.split()[4]))
    if (line.count('Offset')>0) :
      offset[0].append(float(line.split()[2]))
      offset[1].append(float(line.split()[3]))
      offset[2].append(float(line.split()[4]))
    if (line.count('TransformMatrix')>0) :
      axisX[0].append(float(line.split()[2]))
      axisX[1].append(float(line.split()[3]))
      axisX[2].append(float(line.split()[4]))
      axisY[0].append(float(line.split()[5]))
      axisY[1].append(float(line.split()[6]))
      axisY[2].append(float(line.split()[7]))

# compute normalized average direction cosines
# avgZ is computed as the cross product of avgX and avgY
avgX                = [float(np.mean(axisX[0])),  float(np.mean(axisX[1])),  float(np.mean(axisX[2]))]
avgY                = [float(np.mean(axisY[0])),  float(np.mean(axisY[1])),  float(np.mean(axisY[2]))]
avgX                = [avgX[i]/math.sqrt(sum(avgX[i]*avgX[i] for i in range(len(avgX))))  for i in range(len(avgX))]
avgY                = [avgY[i]/math.sqrt(sum(avgY[i]*avgY[i] for i in range(len(avgY))))  for i in range(len(avgY))]
avgZ                = [avgX[1]*avgY[2] - avgX[2]*avgY[1],avgX[2]*avgY[0] - avgX[0]*avgY[2],avgX[0]*avgY[1] - avgX[1]*avgY[0]]

outputDimensions    = [int(np.min(dims[0])), int(np.min(dims[1])), int(np.min(dims[2]))]
outputSpacing       = [float(np.median(spacing[0])), float(np.median(spacing[1])), float(np.median(spacing[2]))]
outputOffset        = [float(np.mean(offset[0])),  float(np.mean(offset[1])),  float(np.mean(offset[2]))]
outputTransform     = [avgX,avgY,avgZ]

# -------------------------------------------------------
# Compute the forward registrations
#
# i.e., each reference dataset is the current fixed image
# -------------------------------------------------------
print("> Running forward registrations...\n")

# check if parameter files exist
if not os.path.isfile(parameterFile0) :
  print (">> Error: parameter file not found...")
  raise SystemExit

if not os.path.isfile(parameterFile1) :
  print (">> Error: parameter file not found...")
  raise SystemExit

# run Elastix
forwardRegistrationSuccess = []

for reference in range(1,dataSetCount+1):
  for subject in range(1,dataSetCount+1):
    if (reference < subject) :
      # create the directory
      registerDir   = "%s_f%s_m%s_%s" % (scalarFile, str(reference).rjust(3, '0'),str(subject).rjust(3, '0'), "fwd")
      resultDir     = os.path.join(outputDir, registerDir)

      if not os.path.exists(resultDir):
        os.makedirs(resultDir)

      # define volumes for current registration
      fixedIndex    = str(reference).rjust(3, '0')
      movingIndex   = str(subject).rjust(3, '0')
      fixedVolume   = os.path.join(currentDir, "%s/%s_%s.mhd" % (dataDir,scalarFile,fixedIndex))
      movingVolume  = os.path.join(currentDir, "%s/%s_%s.mhd" % (dataDir,scalarFile,movingIndex))

      # run Elastix in a subprocess (hiding the large output)
      currentLabel = ">> registration fwd > fixed: " + str(reference).rjust(3, '0') + " - moving: " + str(subject).rjust(3, '0')
      sys.stdout.write(currentLabel)

      elastixCommand     = "elastix -f %s -m %s -out %s -p %s -p %s" % (fixedVolume, movingVolume, resultDir, parameterFile0, parameterFile1)

      if not os.path.exists( os.path.join( resultDir, resultLabel1 )  ) :
        p = subprocess.Popen(elastixCommand.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        p.wait()

        out = out.splitlines()
        success = False
        for line in out:
          if("Total time elapsed" in line):
            seconds = line.replace("Total time elapsed: ", "")
            success = True
      else :
        seconds  = "x"
        success  = True

      if success :
        forwardRegistrationSuccess.append(True)
        print(("\r" + currentLabel + " - success in " + seconds))
      else :
        forwardRegistrationSuccess.append(False)
        print(("\r" + currentLabel + " - failed"))
        break

# --------------------------------------------------------
# Compute inverse registrations
#
# i.e., each reference dataset is the current moving image
# --------------------------------------------------------
print("")
print("> Running inverse registrations...\n")

# check parameter files
if not os.path.isfile(parameterFileInv) :
  print (">> Error: parameter file not found...")
  raise SystemExit

# run Elastix
inverseRegistrationSuccess = []

for reference in range(1,dataSetCount+1):
  for subject in range(1,dataSetCount+1):
    if (reference > subject) :
      # get the forward transform file
      registerDir       = "%s_f%s_m%s_%s" % (scalarFile, str(subject).rjust(3, '0'),str(reference).rjust(3, '0'), "fwd")
      resultDir         = os.path.join(outputDir, registerDir)
      transformFile0    = os.path.join(resultDir, transformLabel0)
      transformFile1    = os.path.join(resultDir, transformLabel1)

      # create the directory
      registerDirInv   = "%s_f%s_m%s_%s" % (scalarFile, str(reference).rjust(3, '0'),str(subject).rjust(3, '0'), "inv")
      resultDirInv     = os.path.join(outputDir, registerDirInv)
      transformFileInv  = os.path.join(resultDirInv, transformLabel0)

      if not os.path.exists(resultDirInv):
        os.makedirs(resultDirInv)

      # define volumes for current registration
      # note that the original fixed volume index (from forward registration) is used
      fixedIndex    = str(subject).rjust(3, '0')
      fixedVolume   = os.path.join(currentDir, "%s/%s_%s.mhd" % (dataDir,scalarFile,fixedIndex))

      # run Elastix in a subprocess (hiding the large output)
      currentLabel = ">> registration inv > fixed: " + str(reference).rjust(3, '0') + " - moving: " + str(subject).rjust(3, '0')
      sys.stdout.write(currentLabel)

      elastixCommand = "elastix -f %s -m %s -out %s -t0 %s -p %s" % (fixedVolume, fixedVolume, resultDirInv, transformFile1, parameterFileInv)

      if not os.path.exists( os.path.join( resultDirInv, transformLabel0 )  ) :
        p = subprocess.Popen(elastixCommand.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        p.wait()

        # make sure no other transform is performed
        if os.path.exists(transformFileInv):
          for line in fileinput.input(transformFileInv, inplace=1):
            print((re.sub(r'\(InitialTransformParametersFileName.+',r'(InitialTransformParametersFileName "NoInitialTransform")', line)), end=' ')

        out = out.splitlines()
        success = False
        for line in out:
          if("Total time elapsed" in line):
            seconds = line.replace("Total time elapsed: ", "")
            success = True
      else :
        seconds  = "x"
        success  = True

      if success :
        inverseRegistrationSuccess.append(True)
        print(("\r" + currentLabel + " - success in " + seconds))
      else :
        inverseRegistrationSuccess.append(False)
        print(("\r" + currentLabel + " - failed"))
        break

# ------------------------------------------------------
# Create the average transformation definitions
#
# i.e., for each reference we create an average
#       transformation to all other subject
# ------------------------------------------------------
print("")
print("> Create average transformations...\n")

forwardAverageSuccess = []

# build a new transformation parameter files for each subject
for reference in range(1,dataSetCount+1) :
  # create the output dir
  registerDirAvg    = "%s%s" % (scalarFile, "_pop_avg")
  resultDirAvg      = os.path.join(outputDir, registerDirAvg)
  transformLabelAvg = "TransformParameters.avg%s.txt" % (str(reference).rjust(3, '0'))
  transformFileIdt  = os.path.join(resultDirAvg, transformLabelIdt)
  transformFileAvg  = os.path.join(resultDirAvg, transformLabelAvg)

  if not os.path.exists(resultDirAvg):
    os.makedirs(resultDirAvg)

  # define volumes for current registration
  fixedIndex    = str(reference).rjust(3, '0')
  fixedVolume   = os.path.join(currentDir, "%s/mhd/%s_%s.mhd" % (fixedIndex,scalarFile,fixedIndex))

  # determine the averaging weights
  weights = ""
  for i in range(dataSetCount):
    weights += str('1 ')

  # determine the subtransforms to all other references
  # i.e., take all transforms where the current reference is the fixed image
  transforms = transformFileIdt
  for subject in range(1,dataSetCount+1):
    transformName = "%s_f%s_m%s" % (scalarFile, str(reference).rjust(3, '0'),str(subject).rjust(3, '0'))
    transformDir  = ""
    if (reference > subject) :
      transformDir = os.path.join(outputDir, transformName) + "_inv"
      transforms = transforms + '\"' + os.path.join(transformDir,transformLabel0) + '\" '
    elif (reference < subject) :
      transformDir = os.path.join(outputDir, transformName) + "_fwd"
      transforms = transforms + '\"' + os.path.join(transformDir,transformLabel1) + '\" '

  # build the identity transform
  createTransformIdentity(transformFileIdt, outputDimensions, outputSpacing, outputOffset, outputTransform)

  # build the average transform from the current reference
  createTransformAverage(transformFileAvg, transforms, dataSetCount, weights, outputDimensions, outputSpacing, outputOffset, outputTransform)

  # check the average transformation on the sphere centers (hide the output)
  # the list of centers is processed by Transformix
  currentLabel = ">> transform centers> reference: " + str(reference).rjust(3, '0')
  sys.stdout.write(currentLabel)

  centersList       = os.path.join(currentDir, "%s/%s.txt"% (dataDir,centersFile))
  transformCommand  = "transformix -def %s -tp %s -out %s -def all" % (centersList, transformFileAvg, resultDirAvg)

  p = subprocess.Popen(transformCommand.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  out, err = p.communicate()
  p.wait()

  out = out.splitlines()
  success = False
  for line in out:
    if("Elapsed time" in line):
      seconds = line.replace("Elapsed time: ", "")
      success = True

  if success :
    forwardAverageSuccess.append(True)
    print(("\r" + currentLabel + " - success in " + seconds))
  else :
    forwardAverageSuccess.append(False)
    print(("\r" + currentLabel + " - failed"))
    break

# ------------------------------------------------------
# Compute inverse average transformation
#
# i.e., average all forward transformations
# (where the reference is the fixed images)
# and subsequently compute the average
# ------------------------------------------------------
print("")
print("> Compute the inverse average transformation...\n")

# check parameter files
if not os.path.isfile(parameterFileInv) :
  print (">> Error: parameter file not found...")
  raise SystemExit

# run Elastix
inverseAverageSuccess = []

for reference in range(1,dataSetCount+1) :
  # get the forward average transform file
  registerDirAvg    = "%s%s" % (scalarFile, "_pop_avg")
  resultDirAvg      = os.path.join(outputDir, registerDirAvg)
  transformLabelAvg = "TransformParameters.avg%s.txt" % (str(reference).rjust(3, '0'))
  transformFileAvg  = os.path.join(resultDirAvg, transformLabelAvg)

  # create the directory
  registerLabelInv = "%s%s" % ('inv', str(reference).rjust(3, '0'))
  registerDirInv   = os.path.join(registerDirAvg, registerLabelInv)
  resultDirInv     = os.path.join(outputDir, registerDirInv)

  if not os.path.exists(resultDirInv):
    os.makedirs(resultDirInv)

  # define volumes for current registration
  # i.e., take all transforms where the current reference is the fixed image
  fixedIndex    = str(reference).rjust(3, '0')
  fixedVolume   = os.path.join(currentDir, "%s/%s_%s.mhd" % (dataDir,scalarFile,fixedIndex))

  # run Elastix in a subprocess (hiding the large output)
  currentLabel = ">> registration inv > reference: " + str(reference).rjust(3, '0')
  sys.stdout.write(currentLabel)

  elastixCommand = "elastix -f %s -m %s -out %s -t0 %s -p %s" % (fixedVolume, fixedVolume, resultDirInv, transformFileAvg, parameterFileInv)

  if not os.path.exists( os.path.join( resultDirInv, transformLabel0 )  ) :
    p = subprocess.Popen(elastixCommand.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    p.wait()

    out = out.splitlines()
    success = False
    for line in out:
      if("Total time elapsed" in line):
        seconds = line.replace("Total time elapsed: ", "")
        success = True
  else :
    seconds  = "x"
    success  = True

  if success :
    inverseAverageSuccess.append(True)
    print(("\r" + currentLabel + " - success in " + seconds))
  else :
    inverseAverageSuccess.append(False)
    print(("\r" + currentLabel + " - failed"))
    break

# run Transformix
transformAverageInverseSuccess = []
transformAverageInverseSeconds = 0

for reference in range(1,dataSetCount+1):
  # get inverse transform file
  registerDirAvg       = "%s%s" % (scalarFile, "_pop_avg")
  resultDirAvg         = os.path.join(outputDir, registerDirAvg)
  registerLabelInv     = "%s%s" % ('inv', str(reference).rjust(3, '0'))
  resultDirInv         = os.path.join(resultDirAvg, registerLabelInv)
  transformFileInv     = os.path.join(resultDirInv, transformLabel0)
  resultDirSubject     = os.path.join(resultDirAvg, str(reference).rjust(3,'0'))

  # define volumes for current registration
  fixedIndex    = str(reference).rjust(3, '0')
  fixedVolume   = os.path.join(currentDir, "%s/%s_%s.mhd" % (dataDir,scalarFile,fixedIndex))

  # create the result dir
  if not os.path.exists(resultDirSubject):
    os.makedirs(resultDirSubject)

  # make sure the correct transformation output data characteristics are set
  if os.path.exists(transformFileInv):
    for line in fileinput.input(transformFileInv, inplace=1):
      print((re.sub(r'\(InitialTransformParametersFileName.+',r'(InitialTransformParametersFileName "NoInitialTransform")', line)), end=' ')
    for line in fileinput.input(transformFileInv, inplace=1):
      print((re.sub(r'\(Size.+',r'(Size '+str(outputDimensions[0])+' '+str(outputDimensions[1])+' '+str(outputDimensions[2])+')', line)), end=' ')
    for line in fileinput.input(transformFileInv, inplace=1):
      print((re.sub(r'\(Spacing.+',r'(Spacing '+str(outputSpacing[0])+' '+str(outputSpacing[1])+' '+str(outputSpacing[2])+')', line)), end=' ')
    for line in fileinput.input(transformFileInv, inplace=1):
      print((re.sub(r'\(Origin.+',r'(Origin '+str(outputOffset[0])+' '+str(outputOffset[1])+' '+str(outputOffset[2])+')', line)), end=' ')
    for line in fileinput.input(transformFileInv, inplace=1):
      print((re.sub(r'\(Direction.+',r'(Direction '+str(outputTransform[0][0])+' '+str(outputTransform[0][1])+' '+str(outputTransform[0][2])+' '+str(outputTransform[1][0])+' '+str(outputTransform[1][1])+' '+str(outputTransform[1][2])+' '+str(outputTransform[2][0])+' '+str(outputTransform[2][1])+' '+str(outputTransform[2][2])+')', line)), end=' ')

  # run Transformix in a subprocess (hiding the large output)
  currentLabel = ">> transformation   > reference: " + str(reference).rjust(3, '0')
  sys.stdout.write(currentLabel)

  transformCommand  = "transformix -in %s -tp %s -out %s -def all" % (fixedVolume, transformFileInv, resultDirSubject)

  p = subprocess.Popen(transformCommand.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  out, err = p.communicate()
  p.wait()

  success = False
  out = out.splitlines()
  for line in out:
    if("Elapsed time" in line):
      seconds = line.replace("Elapsed time: ", "")
      seconds = seconds.replace(" Seconds.", "")
      transformAverageInverseSeconds = int(seconds)
      success = True

  if success :
    transformAverageInverseSuccess.append(True)
    print(("\r" + currentLabel + " - success in " + str(transformAverageInverseSeconds) + " Seconds."))
  else :
    transformAverageInverseSuccess.append(False)
    print(("\r" + currentLabel + " - failed"))

# ------------------------------------------------------
# Compute the average reference subject
# ------------------------------------------------------
print("")
print("> Compute reference subject...\n")

images     = []
dimensions = [1,1,1]
origin     = [0,0,0]

registerDirAvg   = "%s%s" % (scalarFile, "_pop_avg")
resultDirAvg     = os.path.join(outputDir, registerDirAvg)

if(sitkLoaded) :
  for subject in range(1,dataSetCount + 1) :
    resultDirSubject = os.path.join(resultDirAvg, str(subject).rjust(3,'0'))

    reader = sitk.ImageFileReader()
    reader.SetFileName ( os.path.join(resultDirSubject, resultLabel ) )

    volume     = reader.Execute()
    volumeFlat = sitk.GetArrayFromImage( volume ).flatten()

    images.append( volumeFlat )

    if (subject == 1) :
      dimensions = volume.GetSize()
      origin     = volume.GetOrigin()

  referenceFlat = np.mean(images, axis=0)
  reference     = sitk.GetImageFromArray( referenceFlat.reshape(dimensions[2], dimensions[1], dimensions[0]) )
  reference.SetOrigin(origin)

  writer = sitk.ImageFileWriter()
  writer.SetFileName ( os.path.join(resultDirAvg, resultLabel) )
  writer.Execute ( reference )

  print("\r>> voxel average    > success")
else :
  print(("\r" + currentLabel + " - failed: SimpleITK not installed"))