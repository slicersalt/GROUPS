import os, sys
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import csv
import platform
import time
import urllib
import shutil

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
    self.parent.contributors = ["Mahmoud Mostapha (UNC)"] 
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

    # Global variables of the Interface
    # Directories
    self.CollapsibleButton_Directories = self.getWidget('CollapsibleButton_Directories')
    self.RigidAlignmentInputModelsDirectory = self.getWidget('DirectoryButton_RigidAlignmentInputModelsDirectory')
    self.RigidAlignmentInputFiducialFilesDirectory = self.getWidget('DirectoryButton_RigidAlignmentInputFiducialFilesDirectory')
    self.RigidAlignmentCommonSphereDirectory = self.getWidget('DirectoryButton_RigidAlignmentCommonSphereDirectory')
    self.RigidAlignmentOutputSphericalModelsDirectory = self.getWidget('DirectoryButton_RigidAlignmentOutputSphericalModelsDirectory')
    self.RigidAlignmentOutputModelsDirectory = self.getWidget('DirectoryButton_RigidAlignmentOutputModelsDirectory')
    #   Apply CLIs
    self.ApplyButton = self.getWidget('pushButton_RigidAlignment')

    # Connections
    # Directories
    self.RigidAlignmentInputModelsDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    self.RigidAlignmentInputFiducialFilesDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    self.RigidAlignmentCommonSphereDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    self.RigidAlignmentOutputSphericalModelsDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    self.RigidAlignmentOutputModelsDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    #   Apply CLIs
    self.ApplyButton.connect('clicked(bool)', self.onApplyButton)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  # Functions to recover the widget in the .ui file
  def getWidget(self, objectName):
    return self.findWidget(self.widget, objectName)

  def findWidget(self, widget, objectName):
    if widget.objectName == objectName:
      return widget
    else:
      for w in widget.children():
        resulting_widget = self.findWidget(w, objectName)
        if resulting_widget:
          return resulting_widget
    return None
  #
  #   Directories
  #
  def onSelect(self):
    InputModelsDirectory = self.RigidAlignmentInputModelsDirectory.directory.encode('utf-8')
    self.InputModelsDirectory = InputModelsDirectory
    InputFiducialFilesDirectory = self.RigidAlignmentInputFiducialFilesDirectory.directory.encode('utf-8')
    self.InputFiducialFilesDirectory = InputFiducialFilesDirectory
    CommonSphereDirectory = self.RigidAlignmentCommonSphereDirectory.directory.encode('utf-8')
    self.CommonSphereDirectory = CommonSphereDirectory
    OutputSphericalModelsDirectory = self.RigidAlignmentOutputSphericalModelsDirectory.directory.encode('utf-8')
    self.OutputSphericalModelsDirectory = OutputSphericalModelsDirectory
    OutputModelsDirectory = self.RigidAlignmentOutputModelsDirectory.directory.encode('utf-8')
    self.OutputModelsDirectory = OutputModelsDirectory
    # Check if each directory has been choosen
    self.ApplyButton.enabled = self.InputModelsDirectory != "." and self.InputFiducialFilesDirectory != "." and self.CommonSphereDirectory!= "." and self.OutputSphericalModelsDirectory != "." and self.OutputModelsDirectory != "."

  def onApplyButton(self):
    logic = RigidAlignmentModuleLogic()
    endRigidAlignment = logic.runRigidAlignment(modelsDir=self.InputModelsDirectory, fiducialDir=self.InputFiducialFilesDirectory, sphereDir=self.CommonSphereDirectory, outputsphereDir=self.OutputSphericalModelsDirectory, outputsurfaceDir=self.OutputModelsDirectory)
    ## RigidAlignment didn't run because of invalid inputs
    if not endRigidAlignment:
        self.errorLabel.show()

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

  def runRigidAlignment(self, modelsDir, fiducialDir, sphereDir, outputsphereDir, outputsurfaceDir):
      
      # ------------------------------------ # 
      # ---------- RigidAlignment ---------- # 
      # ------------------------------------ # 

      print("--- function runRigidAlignment() ---")
      """
      Calling RigidAlignment CLI
          Arguments:
           --mesh [<std::string> input models directory]
           --landmark [<std::string> input fiducial files directory]
           --sphere [<std::string> common unit sphere]
           --output [<std::string> output sphers directory] 
      """
      print("--- Inspecting Input Data---")
      # List all the vtk files in the modelsDir
      listMesh = os.listdir(modelsDir)
      if listMesh.count(".DS_Store"):
        listMesh.remove(".DS_Store")
      listMesh.sort()      

      # Creation of a file name for the common unit sphere
      listUnitSphere = os.listdir(sphereDir)
      if listUnitSphere.count(".DS_Store"):
        listUnitSphere.remove(".DS_Store")
      #UnitSphere = os.path.join(sphereDir, listUnitSphere[0])
      UnitSphere = sphereDir + '/' + listUnitSphere[0]
      
      print("--- Rigid Alignment Running ---")
      # Creation of the parameters of Rigid Alignment
      RigidAlignment_parameters = {}
      RigidAlignment_parameters["mesh"]       = modelsDir
      RigidAlignment_parameters["landmark"]   = fiducialDir
      RigidAlignment_parameters["sphere"]     = UnitSphere
      RigidAlignment_parameters["output"]     = outputsphereDir
      RA = slicer.modules.rigidalignment
      # Launch Rigid Alignment
      slicer.cli.run(RA, None, RigidAlignment_parameters, wait_for_completion=True)
      print("--- Rigid Alignment Done ---")

      # ------------------------------------ # 
      # ------------ SurfRemesh ------------ # 
      # ------------------------------------ # 
      print("--- function runSurfRemesh() ---")
      """
      Calling SurfRemesh CLI
          Arguments:
           --tempModel [<std::string> input sphere]
           --input [<std::string> input surface]
           --ref [<std::string> common unit sphere]
           --output [<std::string> output surface] 
      """
     
      listSphere = os.listdir(outputsphereDir)
      if listSphere.count(".DS_Store"):
        listSphere.remove(".DS_Store")
      listSphere.sort()

      outputMeshes = []

      for i in range(0,len(listMesh)):
        #Mesh = os.path.join( modelsDir, listMesh[i])
        Mesh = modelsDir + '/' + listMesh[i]
        #Sphere = os.path.join(outputsphereDir, listSphere[i])
        Sphere = outputsphereDir + '/' + listSphere[i]
        #Mesh = os.path.join(outputsurfaceDir, listSphere[i].split("_rotSphere.vtk",1)[0] + '_aligned.vtk')
        OutputMesh = outputsurfaceDir + '/' + listSphere[i].split("_rotSphere.vtk",1)[0] + '_aligned.vtk'
        # Creation of the parameters of SurfRemesh
        SurfRemesh_parameters = {}
        SurfRemesh_parameters["temp"]       = Sphere
        SurfRemesh_parameters["input"]      = Mesh
        SurfRemesh_parameters["ref"]        = UnitSphere
        SurfRemesh_parameters["output"]     = OutputMesh
        SR = slicer.modules.sremesh
        # Launch SurfRemesh
        slicer.cli.run(SR, None, SurfRemesh_parameters, wait_for_completion=True)
        print("--- Surface " + str(i) + " Remesh Done ---")
        # ------------------------------------ # 
        # ------------ Color Maps ------------ # 
        # ------- and encode landmarks ------- #
        # ------------------------------------ # 
        reader_in = vtk.vtkPolyDataReader()
        reader_in.SetFileName(str(Mesh))
        reader_in.Update()
        init_mesh = reader_in.GetOutput()
        phiArray = init_mesh.GetPointData().GetScalars("_paraPhi")

        reader_out = vtk.vtkPolyDataReader()
        reader_out.SetFileName(str(OutputMesh))
        reader_out.Update()
        new_mesh = reader_out.GetOutput()
        new_mesh.GetPointData().SetActiveScalars("_paraPhi")
        new_mesh.GetPointData().SetScalars(phiArray)
        new_mesh.Modified()

        # encode landmarks
        Fiducial = fiducialDir + '/' + listMesh[i].split("_pp_surf_SPHARM",1)[0] + "_fid.fcsv"
        fid = open(Fiducial)
        lines = fid.readlines()
        pts = []
        for line in lines:
          if line[0] != '#':
            s = line.split(',')
            pt = [ float(s[1]), float(s[2]), float(s[3])]
            pts.append(pt)

        loc = vtk.vtkKdTreePointLocator()
        loc.SetDataSet(new_mesh)
        loc.BuildLocator()

        ptArray = vtk.vtkDoubleArray()
        ptArray.SetNumberOfComponents(1)
        ptArray.SetNumberOfValues(new_mesh.GetNumberOfPoints())
        ptArray.SetName('Landmarks')
        for ind in range(0,ptArray.GetNumberOfValues()):
          ptArray.SetValue(ind,0.0)
        
        for l_ind in range(0, len(pts)):
          ind = loc.FindClosestPoint(pts[l_ind])
          ptArray.SetValue(ind,l_ind+1)
        
        new_mesh.GetPointData().AddArray(ptArray)
        # write circle out
        polyDataWriter = vtk.vtkPolyDataWriter()
        polyDataWriter.SetInputData(new_mesh)
        polyDataWriter.SetFileName(str(OutputMesh))
        polyDataWriter.Write()

        outputMeshes.append((new_mesh, str(OutputMesh)))
      print("--- Surf Remesh Done ---")

      print("--- Inspecting Results ---")
      # Load vtk files in ShapePopulationViewer
      slicer.modules.shapepopulationviewer.widgetRepresentation().deleteAll()
      for mesh in outputMeshes:
        polydata, modelName = mesh
        slicer.modules.shapepopulationviewer.widgetRepresentation().loadModel(mesh, modelName)

      slicer.util.selectModule(slicer.modules.shapepopulationviewer)
