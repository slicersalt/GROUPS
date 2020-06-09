import csv
import datetime
import glob
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
    self.parent.title = "SPHARM-PDM Correspondence Improvement"
    self.parent.categories = ["Shape Creation"]
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
    self.ui.CommonSphereDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    self.ui.FiducialsDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    self.ui.OutputDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    self.ui.OutputSphereDirectory.connect('directoryChanged(const QString &)', self.onSelect)
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
    self.inputDir = pathlib.Path(self.ui.InputDirectory.directory)
    self.commonSphereDir = pathlib.Path(self.ui.CommonSphereDirectory.directory)
    self.fiducialsDir = pathlib.Path(self.ui.FiducialsDirectory.directory)
    self.outputDir = pathlib.Path(self.ui.OutputDirectory.directory)
    self.outputSphereDir = pathlib.Path(self.ui.OutputSphereDirectory.directory)

    # Check if each directory has been choosen
    self.ui.ApplyButton.enabled = '.' not in (self.inputDir, self.fiducialsDir, self.outputDir)

  def onApplyButton(self):
    models = self.inputDir.glob('*_pp_surf_SPHARM.vtk')
    fiducials = self.fiducialsDir.glob('*_fid.fcsv')
    unitSphere = next(self.commonSphereDir.glob('*_surf_para.vtk'))

    logic = RigidAlignmentModuleLogic()
    logic.run(
      models=models,
      fiducials=fiducials,
      unitSphere=unitSphere,
      outModelsDir=self.outputDir,
      outSphereDir=self.outputSphereDir,
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

  def run(self, models, fiducials, unitSphere, outModelsDir, outSphereDir):
    """
    models: A sequence of paths to SPHARM model files. (*_pp_surf_SPHARM.vtk)
    fiducials: A sequence of paths to fiducial data files. (*_fid.fcsv)
    unitSphere: A path to a unit sphere for alignment. (*_surf_para.vtk)
    outputDir: Output directory for aligned spheres.
    """

    models = sorted(models)
    fiducials = sorted(fiducials)
    body = list(zip(models, fiducials))

    temp = pathlib.Path(slicer.util.tempDirectory(key='RigidAlignment'))

    now = datetime.datetime.now().isoformat()
    inputCSV = temp / '{}.csv'.format(now)

    with inputCSV.open('w', newline='') as f:
      for row in body:
        row = (str(e) + ',' for e in row)
        line = ''.join(row) + '\n'
        f.write(line)

    self.runRigidAlignment(inputCSV, unitSphere, outSphereDir)

    results = []

    for model, fiducial in body:
      name = model.name.rsplit('_pp_surf_SPHARM', 1)[0]

      sphere = os.path.join(outSphereDir, name + '_rotSphere.vtk')
      outModel = os.path.join(outModelsDir, name + '_aligned.vtk')

      self.runSurfRemesh(sphere, model, unitSphere, outModel)
      res = self.buildColorMap(model, fiducial, outModel)

      results.append(res)

    if results:
      self.showViewer(results)

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
      'temp': str(sphere),
      'input': str(model),
      'ref': str(unitSphere),
      'output': str(outModel),
    }

    logging.info('Launching SRemesh Module.')
    slicer.cli.run(slicer.modules.sremesh, None, args, wait_for_completion=True)
    logging.info('SRemesh Completed.')

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
