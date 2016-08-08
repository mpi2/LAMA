surfImg = vtk.vtkPolyDataReader()
surfImg.SetFileName('/home/henrik/Downloads/test_skeleton_coords.vtk')

pd = vtk.vtkPolyData()
pd = surfImg.GetOutput()
pd.Update()

# spheres

sphere = vtk.vtkSphereSource()
sphere.SetRadius(0.25)
sphere.SetThetaResolution(25)
sphere.SetPhiResolution(25)

glyphs = vtk.vtkGlyph3D()
glyphs.SetInput(pd)
glyphs.SetSource(sphere.GetOutput())
glyphs.ScalingOff()

# tubes

numOfPoints = pd.GetNumberOfPoints()

lines = vtk.vtkCellArray()
linepoints = vtk.vtkPoints()

j = 0
for i in range(0, numOfPoints, 1):
    p = pd.GetPoint(i)
    p2 = pd.GetPoint(i + 1)
lines.InsertNextCell(2)
lines.InsertCellPoint(j)
linepoints.InsertPoint(j, p[0], p[1], p[2])
j = j + 1
lines.InsertCellPoint(j)
linepoints.InsertPoint(j, p2[0], p2[1], p2[2])

tubeData = vtk.vtkPolyData()
tubeData.SetPoints(linepoints)
tubeData.SetLines(lines)

Tubes = vtk.vtkTubeFilter()
Tubes.SetInput(tubeData)
Tubes.SetRadius(0.1)
Tubes.SetNumberOfSides(10)

model = slicer.vtkMRMLModelNode()
model.SetAndObservePolyData(glyphs.GetOutput())
modelDisplay = slicer.vtkMRMLModelDisplayNode()
slicer.mrmlScene.AddNode(modelDisplay)
model.SetAndObserveDisplayNodeID(modelDisplay.GetID())
modelDisplay.SetInputPolyData(model.GetPolyData())
slicer.mrmlScene.AddNode(model)

model = slicer.vtkMRMLModelNode()
model.SetAndObservePolyData(Tubes.GetOutput())
modelDisplay = slicer.vtkMRMLModelDisplayNode()
slicer.mrmlScene.AddNode(modelDisplay)
model.SetAndObserveDisplayNodeID(modelDisplay.GetID())
modelDisplay.SetInputPolyData(model.GetPolyData())
slicer.mrmlScene.AddNode(model)
