import vtk



"""
Scales an organ mess and spits out vtk file
"""
SCALE = 0.5


filename = "/home/neil/work/fake-abnormalities/out/1.vtk"

reader = vtk.vtkPolyDataReader()
reader.SetFileName(filename)
reader.Update()

poly = reader.GetOutput()

com_filter = vtk.vtkCenterOfMass()
com_filter.SetInputData(poly)
com_filter.SetUseScalarsAsWeights(False)
com_filter.Update()

cx, cy, cz = com_filter.GetCenter()


#Do transform

tform = vtk.vtkTransform()
tform.Translate(cx, cy, cz)
tform.Scale(SCALE, SCALE, SCALE)
tform.Translate(-cx, -cy, -cz)

tfilter = vtk.vtkTransformPolyDataFilter()
tfilter.SetInputData(poly)
tfilter.SetTransform(tform)
tfilter.Update()

out = 'scaled_mesh.vtk'
writer = vtk.vtkPolyDataWriter()
writer.SetInputConnection(tfilter.GetOutputPort())
writer.SetFileName(out)
writer.Write()
