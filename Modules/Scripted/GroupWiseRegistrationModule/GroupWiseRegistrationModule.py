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
# GroupWiseRegistrationModule
#

class GroupWiseRegistrationModule(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
      ScriptedLoadableModule.__init__(self, parent)
      self.parent.title = "Group Wise Registeration Module"
      self.parent.categories = ["Groups"]
      self.parent.dependencies = []
      self.parent.contributors = ["Mahmoud Mostapha (UNC)"]
      self.parent.helpText = """
      Cortical correspondence method employing group-wise registration in a spherical parametrization space for the use in neuroimaging studies.
      The proposed method is unbiased registration that estimates a continuous smooth deformation field into an unbiased average space via entropy minimization
      using spherical harmonic decomposition of the spherical deformation field.
      """
      self.parent.acknowledgementText = """
        This work was supported by NIH NIBIB R01EB021391
        (Shape Analysis Toolbox for Medical Image Computing Projects).
      """

#
# GroupWiseRegistrationModuleWidget
#

class GroupWiseRegistrationModuleWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)
    #
    #  Interface
    #
    loader = qt.QUiLoader()
    self.moduleName = 'GroupWiseRegistrationModule'
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
    self.GroupsInputModelsDirectory = self.getWidget('DirectoryButton_GroupsInputModelsDirectory')
    self.GroupsInputSphericalModelsDirectory = self.getWidget('DirectoryButton_GroupsInputSphericalModelsDirectory')
    self.GroupsOutputCoefficientsDirectory = self.getWidget('DirectoryButton_GroupsOutputCoefficientsDirectory')
    self.GroupsOutputModelsDirectory = self.getWidget('DirectoryButton_GroupsOutputModelsDirectory')
    # Parameters
    self.CollapsibleButton_Parameters = self.getWidget('CollapsibleButton_Parameters')
    self.GroupsLandmarks = self.getWidget('checkBox_Landmarks')
    self.GroupsSpharmDegree = self.getWidget('SliderWidget_SpharmDecomposotion')
    self.GroupsIterations = self.getWidget('spinBox_Iterations')
    self.GroupsProperties = self.getWidget('tableWidget_Properites')
    # Apply CLIs
    self.ApplyButton = self.getWidget('pushButton_RunGroups')

    # Initial Values
    self.Landmarks = False
    self.SpharmDegree = 5
    self.Iterations = 5000
    self.PropertiesNames = []
    self.PropertiesNamesSelected = []
    self.PropertiesWeightsSelected = []

    # Connections
    # Directories
    self.CollapsibleButton_Directories.connect('clicked()',
                                                        lambda: self.onSelectedCollapsibleButtonOpen(
                                                          self.CollapsibleButton_Directories))
    self.GroupsInputModelsDirectory.connect('directoryChanged(const QString &)', self.onSelectInputModels)
    self.GroupsInputSphericalModelsDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    self.GroupsOutputCoefficientsDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    self.GroupsOutputModelsDirectory.connect('directoryChanged(const QString &)', self.onSelect)
    # Parameters
    self.CollapsibleButton_Parameters.connect('clicked()',
                                                        lambda: self.onSelectedCollapsibleButtonOpen(
                                                          self.CollapsibleButton_Parameters))
    self.GroupsLandmarks.connect('stateChanged(int)', self.onSelectLandmarks)
    self.GroupsSpharmDegree.connect('valueChanged(double)', self.onSelectSpharmDegree)
    self.GroupsIterations.connect('valueChanged(int)', self.onSelectIterations)
    # Apply CLIs
    self.ApplyButton.connect('clicked(bool)', self.onApplyButton)

    # Widget Configuration
    # Input Properities Table Configuration
    self.GroupsProperties.setColumnCount(2)
    self.GroupsProperties.setHorizontalHeaderLabels([' Properity ', ' Weight '])
    self.GroupsProperties.setColumnWidth(0, 400)
    horizontalHeader = self.GroupsProperties.horizontalHeader()
    horizontalHeader.setStretchLastSection(False)
    horizontalHeader.setResizeMode(0, qt.QHeaderView.Stretch)
    horizontalHeader.setResizeMode(1, qt.QHeaderView.ResizeToContents)


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

  # Only one tab can be displayed at the same time:
  # When one tab is opened all the other tabs are closed
  def onSelectedCollapsibleButtonOpen(self, selectedCollapsibleButton):
    if selectedCollapsibleButton.isChecked():
      collapsibleButtonList = [self.CollapsibleButton_Directories,
                               self.CollapsibleButton_Parameters]
      for collapsibleButton in collapsibleButtonList:
        collapsibleButton.setChecked(False)
      selectedCollapsibleButton.setChecked(True)

  def onSelectInputModels(self):
    InputModelsDirectory = self.GroupsInputModelsDirectory.directory.encode('utf-8')
    self.InputModelsDirectory = InputModelsDirectory
    self.PropertiesNames = []
    listMesh = os.listdir(InputModelsDirectory)
    if listMesh.count(".DS_Store"):
      listMesh.remove(".DS_Store")
    #VTKfilepath = os.path.join( modelsDir, listMesh[0])
    VTKfilepath = InputModelsDirectory + '/' + listMesh[0]
    if os.path.exists(VTKfilepath):
      reader = vtk.vtkPolyDataReader()
      reader.SetFileName(VTKfilepath)
      reader.Update()
      polydata = reader.GetOutput()
      for i in range(polydata.GetPointData().GetNumberOfArrays()):
        self.PropertiesNames.append(polydata.GetPointData().GetArrayName(i))
    row = 0
    for PropertyName in self.PropertiesNames:
      self.GroupsProperties.setRowCount(row + 1)
      # Column 0:
      labelPropertyName = qt.QLabel(PropertyName)
      labelPropertyName.setAlignment(0x84)
      self.GroupsProperties.setCellWidget(row, 0, labelPropertyName)
      # Column 1:
      widget = qt.QWidget()
      layout = qt.QHBoxLayout(widget)
      spinBox = qt.QSpinBox()
      spinBox.setMinimum(0)
      layout.addWidget(spinBox)
      layout.setAlignment(0x84)
      layout.setContentsMargins(0, 0, 0, 0)
      widget.setLayout(layout)
      self.GroupsProperties.setCellWidget(row, 1, widget)
      row = row + 1

  def onSelect(self):
    InputSphericalModelsDirectory = self.GroupsInputSphericalModelsDirectory.directory.encode('utf-8')
    self.InputSphericalModelsDirectory = InputSphericalModelsDirectory
    OutputCoefficientsDirectory = self.GroupsOutputCoefficientsDirectory.directory.encode('utf-8')
    self.OutputCoefficientsDirectory = OutputCoefficientsDirectory
    OutputModelsDirectory = self.GroupsOutputModelsDirectory.directory.encode('utf-8')
    self.OutputModelsDirectory = OutputModelsDirectory
    # Check if each directory has been choosen
    self.ApplyButton.enabled = self.InputModelsDirectory != "." and self.InputSphericalModelsDirectory != "." and self.OutputCoefficientsDirectory!= "." and self.OutputModelsDirectory != "."

  def onSelectLandmarks(self):
    if self.GroupsLandmarks.checkState():
      self.Landmarks = True
    else:
      self.Landmarks = False

  def onSelectSpharmDegree(self):
    self.SpharmDegree = self.GroupsSpharmDegree.value

  def onSelectIterations(self):
    self.Iterations = self.GroupsIterations.value

  def onApplyButton(self):

    self.PropertiesNamesSelected = []
    self.PropertiesWeightsSelected = []
    # Check Selected Properities ( > 0 )
    table = self.GroupsProperties
    for row in range(table.rowCount):
        widget = table.cellWidget(row, 1)
        tuple = widget.children()
        spinbox = tuple[1]
        if (spinbox.value > 0):
          self.PropertiesNamesSelected.append(self.PropertiesNames[row])
          self.PropertiesWeightsSelected.append(spinbox.value)

    logic = GroupWiseRegistrationModuleLogic()
    endGroups = logic.runGroups(modelsDir=self.InputModelsDirectory, sphereDir=self.InputSphericalModelsDirectory, outputCoefficientDir=self.OutputCoefficientsDirectory, outputDir=self.OutputModelsDirectory, Landmarks=self.Landmarks,
                                    propertiesNames = self.PropertiesNamesSelected, propertiesWeights = self.PropertiesWeightsSelected, degree = self.SpharmDegree, maxIter = self.Iterations)
#
# GroupWiseRegistrationModuleLogic
#
class GroupWiseRegistrationModuleLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def runGroups(self, modelsDir, sphereDir, outputCoefficientDir, outputDir, Landmarks=False, propertiesNames=[], propertiesWeights=[], degree=5, maxIter=5000):
    print "--- function runGroupWiseRegisteration() ---"
    """
    Calling Groups CLI
        Arguments:
         --surface: Directory with input models
         --sphere: Sphere folder
         --output: Output directory
         --modelProperty: Use a property embeded in the VTK file (need fixing)
         --landmarksOn:  Activate the use of landmarks
         --degree: Degree of deformation field
         --maxIter: Maximum number of iteration
    """
    print "--- Inspecting Input Data---"
    # List all the vtk files in the modelsDir
    listMesh = os.listdir(modelsDir)
    if listMesh.count(".DS_Store"):
      listMesh.remove(".DS_Store")
    # Creation of a CSV file to load the vtk files in ShapePopulationViewer
    #filePathCSV = os.path.join( slicer.app.temporaryPath, 'PreviewInputDataForVisualizationInSPV.csv')
    filePathCSV = slicer.app.temporaryPath + '/' + 'PreviewInputDataForVisualizationInSPV.csv'
    file = open(filePathCSV, 'w')
    cw = csv.writer(file, delimiter=',')
    cw.writerow(['VTK Files'])
    # Add the path of the vtk files
    for i in range(0, len(listMesh)):
      #VTKfilepath = os.path.join( modelsDir, listMesh[i])
      VTKfilepath = modelsDir + '/' + listMesh[i]
      if os.path.exists(VTKfilepath):
        cw.writerow([VTKfilepath])
    file.close()
    # Creation of the parameters of SPV
    parameters = {}
    parameters["CSVFile"] = filePathCSV
    #   If a binary of SPV has been installed
    if hasattr(slicer.modules, 'shapepopulationviewer'):
      SPV = slicer.modules.shapepopulationviewer
    #   If SPV has been installed via the Extension Manager
    elif hasattr(slicer.modules, 'launcher'):
      SPV = slicer.modules.launcher
    # Launch SPV
    slicer.cli.run(SPV, None, parameters, wait_for_completion=True)

    # Deletion of the CSV files in the Slicer temporary directory
    if os.path.exists(filePathCSV):
      os.remove(filePathCSV)
    print "--- Groups Running ---"
    # Creation of the parameters of Rigid Alignment
    Groups_parameters = {}
    Groups_parameters["dirSurf"]           = modelsDir
    Groups_parameters["dirSphere"]         = sphereDir
    Groups_parameters["dirOutput"]         = outputCoefficientDir

    Prop = propertiesNames[0] + "," + str(propertiesWeights[0])
    if len(propertiesNames) > 1 :
      i = 1
      while i < len(propertiesNames):
        Prop = Prop + "," + propertiesNames[i] + "," + str(propertiesWeights[i])
        i += 1
    Groups_parameters["modelProperty"]  = Prop
    if (Landmarks):
      Groups_parameters["landmarksOn"]  = " "
    Groups_parameters["degree"]         = degree
    Groups_parameters["maxIter"]        = maxIter
    #print(Groups_parameters)
    G = slicer.modules.groups
    # Launch Groups
    cliNode = slicer.cli.run(G, None, Groups_parameters, wait_for_completion=True)
    def printStatus(caller, event):
      print("Got a %s from a %s" % (event, caller.GetClassName()))
      if caller.IsA('vtkMRMLCommandLineModuleNode'):
        print("Status is %s" % caller.GetStatusString())
    cliNode.AddObserver('ModifiedEvent', printStatus)
    print "--- Groups Done ---"

    # ------------------------------------ #
    # ------------ SurfRemesh ------------ #
    # ------------------------------------ #
    print "--- function runSurfRemesh() ---"
    """
    Calling SurfRemesh CLI
        Arguments:
         --tempModel : Template model (unit sphere) used for the spherical parameterization
         --input     : Subject surface model
         --ref       : Reference model (unit sphere) for the final surface remeshing
         --coeff     : Spherical harmonic coefficients
         --output    : output surface model
    """
    listSphere = os.listdir(sphereDir)
    if listSphere.count(".DS_Store"):
      listSphere.remove(".DS_Store")

    listCoeff = os.listdir(outputCoefficientDir)
    if listCoeff.count(".DS_Store"):
      listCoeff.remove(".DS_Store")

    for i in range(0,len(listMesh)):
      #Mesh = os.path.join( modelsDir, listMesh[i])
      Mesh = modelsDir + '/' + listMesh[i]
      #Sphere = os.path.join(outputsphereDir, listSphere[i])
      Sphere = sphereDir + '/' + listSphere[i]
      #Coeff = os.path.join(outputCoefficientDir, listCoeff[i])
      Coeff = outputCoefficientDir + '/' + listCoeff[i]
      #OutputMesh = os.path.join(outputDir, listCoeff[i].split(".coeff.txt",1)[0] + '.vtk')
      OutputMesh = outputDir + '/' + listCoeff[i].split(".coeff",1)[0] + '.vtk'
      # Creation of the parameters of SurfRemesh
      SurfRemesh_parameters = {}
      SurfRemesh_parameters["temp"]       = Sphere
      SurfRemesh_parameters["input"]      = Mesh
      SurfRemesh_parameters["ref"]        = Sphere
      SurfRemesh_parameters["coeff"]      = Coeff
      SurfRemesh_parameters["output"]     = OutputMesh
      SR = slicer.modules.sremesh
      # Launch SurfRemesh
      slicer.cli.run(SR, None, SurfRemesh_parameters, wait_for_completion=True)
      print "--- Surface " + str(i) + " Remesh Done ---"
      # ------------------------------------ #
      # ------------ Color Maps ------------ #
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
      # write results
      polyDataWriter = vtk.vtkPolyDataWriter()
      polyDataWriter.SetInputData(new_mesh)
      polyDataWriter.SetFileName(str(OutputMesh))
      polyDataWriter.Write()
    print "--- Surf Remesh Done ---"

    print "--- Inspecting Results ---"
    # List all the vtk files in the outputDir
    listOutputMesh = os.listdir(outputDir)
    if listOutputMesh.count(".DS_Store"):
      listOutputMesh.remove(".DS_Store")
    # Creation of a CSV file to load the output vtk files in ShapePopulationViewer
    #filePathCSV = os.path.join( slicer.app.temporaryPath, 'PreviewOutputDataForVisualizationInSPV.csv')
    filePathCSV = slicer.app.temporaryPath + '/' + 'PreviewOutputDataForVisualizationInSPV.csv'
    file = open(filePathCSV, 'w')
    cw = csv.writer(file, delimiter=',')
    cw.writerow(['VTK Files'])
    # Add the path of the vtk files
    for i in range(0, len(listOutputMesh)):
      #VTKfilepath = os.path.join( outputDir, listOutputMesh[i])
      VTKfilepath = outputDir + '/' + listOutputMesh[i]
      if os.path.exists(VTKfilepath):
        cw.writerow([VTKfilepath])
    file.close()

    # Creation of the parameters of SPV
    parameters = {}
    parameters["CSVFile"] = filePathCSV
    # Launch SPV
    slicer.cli.run(SPV, None, parameters, wait_for_completion=True)

    # Deletion of the CSV files in the Slicer temporary directory
    if os.path.exists(filePathCSV):
      os.remove(filePathCSV)
