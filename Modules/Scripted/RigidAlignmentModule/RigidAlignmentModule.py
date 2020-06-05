import os
import sys
import unittest
import uuid

import vtk
import qt
import ctk
import slicer
import slicer.util
from slicer.ScriptedLoadableModule import *
import logging
import csv
import platform
import time
import urllib
import shutil
import glob
import datetime


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
    self.ui.FiducialsDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    self.ui.OutputDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    #   Apply CLIs
    self.ui.ApplyButton.connect('clicked(bool)', self.onApplyButton)

    self.ui.InputDirectory.directory = '/home/allem/Downloads/groups/models'
    self.ui.FiducialsDirectory.directory = '/home/allem/Downloads/groups/fiducials'
    self.ui.OutputDirectory.directory = '/home/allem/Downloads/groups/out'

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

def pluckPattern(pattern):
  """Each file matched by the glob pattern will be mapped to by a hard link in the temporary output directory. The
  directory path is returned. """

  key = str(uuid.uuid4())
  target = slicer.util.tempDirectory(key=key)

  for src in glob.iglob(pattern):
    name = os.path.basename(src)
    dst = os.path.join(target, name)
    os.link(src, dst)

  return target


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
    models = pluckPattern(os.path.join(modelsDir, '*_pp_surf_SPHARM.vtk'))
    fiducials = pluckPattern(os.path.join(fiducialsDir, '*_fid.fcsv'))
    unitSphere = next(glob.iglob(os.path.join(sphereDir, '*_surf_para.vtk')))

    results = []

    try:
      self.runRigidAlignment(models, fiducials, unitSphere, outSphereDir)

      for name in os.listdir(models):
        name = name.rsplit('_pp_surf_SPHARM', 1)[0]

        model = os.path.join(models, name + '_pp_surf_SPHARM.vtk')
        sphere = os.path.join(outSphereDir, name + '_rotSphere.vtk')
        outModel = os.path.join(outModelsDir, name + '_aligned.vtk')
        fiducial = os.path.join(fiducials, name + '_fid.fcsv')

        self.runSurfRemesh(sphere, model, unitSphere, outModel)
        res = self.buildColorMap(model, fiducial, outModel)

        results.append(res)
    finally:
      shutil.rmtree(models)
      shutil.rmtree(fiducials)

    self.showViewer(results)

  def runRigidAlignment(self, models, fiducials, sphere, outputDir):
    args = {
      'mesh': models,
      'landmark': fiducials,
      'sphere': sphere,
      'output': outputDir
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
