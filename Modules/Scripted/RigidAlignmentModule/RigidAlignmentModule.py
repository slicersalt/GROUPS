import csv
import datetime
import logging
import os
import pathlib

import qt
import slicer
import slicer.util
import slicer.cli
from slicer.ScriptedLoadableModule import *
import vtk


#
# RigidAlignmentModule
#

class RigidAlignmentModule(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Rigid Alignment Module" 
    self.parent.categories = ["Groups"]
    self.parent.dependencies = []
    self.parent.contributors = ["Mahmoud Mostapha (UNC), Jared Vicory (Kitware), David Allemang (Kitware)"]
    self.parent.helpText = """
    Rigid alignment of the landmarks on the unit sphere: the input models share the same unit sphere 
    and their landmarks are defined as spacial coordinates (x,y,z) of the input model. 
    """
    self.parent.acknowledgementText = """
      This work was supported by NIH NIBIB R01EB021391
      (Shape Analysis Toolbox for Medical Image Computing Projects).
    """


#
# RigidAlignmentModuleWidget
#

class RigidAlignmentModuleWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)
    #
    #  Interface
    #
    loader = qt.QUiLoader()
    self.moduleName = 'RigidAlignmentModule'
    scriptedModulesPath = eval('slicer.modules.%s.path' % self.moduleName.lower())
    scriptedModulesPath = os.path.dirname(scriptedModulesPath)
    path = os.path.join(scriptedModulesPath, 'Resources', 'UI', '%s.ui' % self.moduleName)
    qfile = qt.QFile(path)
    qfile.open(qt.QFile.ReadOnly)
    widget = loader.load(qfile, self.parent)
    self.layout = self.parent.layout()
    self.widget = widget
    self.layout.addWidget(widget)
    self.ui = slicer.util.childWidgetVariables(widget)

    # Connections
    # Directories
    self.ui.InputDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    self.ui.FiducialsDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    self.ui.OutputDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    #   Apply CLIs
    self.ui.ApplyButton.connect('clicked(bool)', self.onApplyButton)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  #
  #   Directories
  #
  def onSelect(self):
    self.inputDir = self.ui.InputDirectory.directory
    self.fiducialsDir = self.ui.FiducialsDirectory.directory
    self.outputDir = self.ui.OutputDirectory.directory

    # Check if each directory has been choosen
    self.ui.ApplyButton.enabled = '.' not in (self.inputDir, self.fiducialsDir, self.outputDir)

  def onApplyButton(self):
    logic = RigidAlignmentModuleLogic()
    logic.run(
      modelsDir=self.inputDir,
      sphereDir=self.inputDir,
      fiducialsDir=self.fiducialsDir,
      outModelsDir=self.outputDir,
      outSphereDir=self.outputDir,
    )

#
# RigidAlignmentModuleLogic
#

class RigidAlignmentModuleLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def run(self, modelsDir, sphereDir, fiducialsDir, outModelsDir, outSphereDir):
    unitSphere = next(pathlib.Path(sphereDir).glob('*_surf_para.vtk'))
    models = pathlib.Path(modelsDir).glob('*_pp_surf_SPHARM.vtk')
    fiducials = pathlib.Path(fiducialsDir).glob('*_fid.fcsv')

    temp = pathlib.Path(slicer.util.tempDirectory(key='RigidAlignment'))

    now = datetime.datetime.now().isoformat()
    inputCSV = temp / '{}.csv'.format(now)
    body = zip(sorted(fiducials), sorted(models))

    with inputCSV.open('w') as f:
      writer = csv.writer(f)
      writer.writerows(body)

    self.runRigidAlignment(inputCSV, unitSphere, outSphereDir)

  def runRigidAlignment(self, inputCSV, sphere, outputDir):
    args = {
      'inputCSV': str(inputCSV),
      'sphere': str(sphere),
      'output': str(outputDir)
    }

    logging.info('Launching RigidAlignment Module.')
    slicer.cli.run(slicer.modules.rigidalignment, None, args, wait_for_completion=True)
    logging.info('RigidAlignment Completed.')

  def runSurfRemesh(self, sphere, model, unitSphere, outModel):
    args = {
      'temp': sphere,
      'input': model,
      'ref': unitSphere,
      'output': outModel,
    }

    try:
      logging.info('Launching SRemesh Module.')
      slicer.cli.run(slicer.modules.sremesh, None, args, wait_for_completion=True)
      logging.info('SRemesh Completed.')
    finally:
      pass

  def buildColorMap(self, model, fiducial, outModel):
    reader_in = vtk.vtkPolyDataReader()
    reader_in.SetFileName(str(model))
    reader_in.Update()
    init_mesh = reader_in.GetOutput()
    phiArray = init_mesh.GetPointData().GetScalars("_paraPhi")

    reader_out = vtk.vtkPolyDataReader()
    reader_out.SetFileName(str(outModel))
    reader_out.Update()
    new_mesh = reader_out.GetOutput()
    new_mesh.GetPointData().SetActiveScalars("_paraPhi")
    new_mesh.GetPointData().SetScalars(phiArray)
    new_mesh.Modified()

    with open(fiducial) as fid:
      lines = fid.readlines()
      pts = []
      for line in lines:
        if line[0] != '#':
          s = line.split(',')
          pt = [float(s[1]), float(s[2]), float(s[3])]
          pts.append(pt)

    loc = vtk.vtkKdTreePointLocator()
    loc.SetDataSet(new_mesh)
    loc.BuildLocator()

    ptArray = vtk.vtkDoubleArray()
    ptArray.SetNumberOfComponents(1)
    ptArray.SetNumberOfValues(new_mesh.GetNumberOfPoints())
    ptArray.SetName('Landmarks')
    for ind in range(0, ptArray.GetNumberOfValues()):
      ptArray.SetValue(ind, 0.0)

    for l_ind in range(0, len(pts)):
      ind = loc.FindClosestPoint(pts[l_ind])
      ptArray.SetValue(ind, l_ind + 1)

    new_mesh.GetPointData().AddArray(ptArray)

    # write results
    polyDataWriter = vtk.vtkPolyDataWriter()
    polyDataWriter.SetInputData(new_mesh)
    polyDataWriter.SetFileName(str(outModel))
    polyDataWriter.Write()

    return new_mesh, str(outModel)

  def showViewer(self, results):
    viewer = slicer.modules.shapepopulationviewer.widgetRepresentation()
    viewer.deleteModels()
    for polydata, name in results:
      viewer.loadModel(polydata, name)
    slicer.util.selectModule(slicer.modules.shapepopulationviewer)
